version 1.0
import "../../tasks/quality_control/read_filtering/task_porechop.wdl" as porechop_task
import "../../tasks/quality_control/read_filtering/task_nanoq.wdl" as nanoq_task
import "../../tasks/quality_control/basic_statistics/task_nanoplot.wdl" as nanoplot_task
import "../../tasks/quality_control/basic_statistics/task_quast.wdl" as quast_task
import "../../tasks/quality_control/basic_statistics/task_assembly_metrics.wdl" as assembly_metrics_task
import "../../tasks/quality_control/basic_statistics/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/quality_control/advanced_metrics/task_checkv.wdl" as checkv_task
import "../../tasks/quality_control/comparisons/task_screen.wdl" as screen_task
import "../../tasks/taxon_id/contamination/task_metabuli.wdl" as metabuli_task
import "../../tasks/taxon_id/task_skani.wdl" as skani_task
import "../../tasks/utilities/task_rasusa.wdl" as rasusa_task
import "../../tasks/taxon_id/task_ete4_taxon_id.wdl" as identify_taxon_id_task
import "../../tasks/utilities/data_handling/task_parse_mapping.wdl" as parse_mapping_task
import "../../tasks/utilities/data_handling/task_fasta_utilities.wdl" as fasta_utilities_task
import "../../tasks/assembly/task_raven.wdl" as raven_task
import "../../tasks/assembly/task_flye.wdl" as flye_task
import "../../tasks/assembly/task_bcftools_consensus.wdl" as bcftools_consensus_task
import "../../tasks/alignment/task_minimap2.wdl" as minimap2_task
import "../../tasks/gene_typing/variant_detection/task_clair3_variants.wdl" as clair3_task
import "../../tasks/task_versioning.wdl" as versioning_task
import "../../tasks/quality_control/read_filtering/task_ncbi_scrub.wdl" as ncbi_scrub_task
import "../utilities/wf_host_decontaminate.wdl" as host_decontaminate_wf
import "../utilities/wf_morgana_magic.wdl" as morgana_magic_wf

workflow theiaviral_ont {
  meta {
    description: "De novo assembly, dynamic reference selection, and reference-based consensus calling for viral metagenomic/genomic data generated on ONT NGS platforms."
  }
  input {
    File read1
    String taxon #taxon id OR organism name (both work)
    String read_extraction_rank = "family"
    String samplename
    String? host # host genome to dehost reads, if desired
    File? reference_fasta
    File? reference_gene_locations_bed # optional, if provided will be used for coverage calculations
    Int? genome_length
    Int min_depth = 10 # minimum depth for masking low coverage regions
    Float min_allele_freq = 0.6 # minimum allele frequency for consensus calling
    Int min_map_quality = 20 # minimum read mapping quality
    Boolean extract_unclassified = true
    Boolean call_porechop = false
    Boolean skip_rasusa = true 
    Boolean skip_screen = false
    Boolean call_raven = true
  }
  # get the PHB version
  call versioning_task.version_capture {
    input:
  }
  # get the taxon id, taxon name, and taxon rank from the user provided taxon
  call identify_taxon_id_task.ete4_taxon_id as ete4_identify {
    input:
      taxon = taxon,
      rank = read_extraction_rank
  }
  # raw read quality check
  call nanoplot_task.nanoplot as nanoplot_raw {
    input:
      read1 = read1,
      samplename = samplename,
      est_genome_length = select_first([genome_length, 12500]) # default viral genome length
  }
  # adapter trimming
  if (call_porechop) {
    call porechop_task.porechop as porechop {
      input:
        read1 = read1,
        samplename = samplename
    }
  }
  # read filtering
  call nanoq_task.nanoq as nanoq {
    input:
      read1 = select_first([porechop.trimmed_reads, read1]),
      samplename = samplename
  }
  # human read scrubbing
  call ncbi_scrub_task.ncbi_scrub_se {
    input:
      read1 = nanoq.filtered_read1,
      samplename = samplename,
  }
  # decontaminate host reads if a host genome is provided
  if (defined(host)) {
    call host_decontaminate_wf.host_decontaminate {
      input:
        samplename = samplename,
        read1 = ncbi_scrub_se.read1_dehosted,
        host = select_first([host])
    }
  }
  # taxonomic classification and read extraction
  call metabuli_task.metabuli as metabuli {
    input:
      read1 = select_first([host_decontaminate.dehost_read1, ncbi_scrub_se.read1_dehosted]),
      samplename = samplename,
      taxon_id = select_first([ete4_identify.taxon_id]),
      extract_unclassified = extract_unclassified
  }
  # downsample reads if the user wants, rasusa parameters are set in the task
  if (! skip_rasusa) {
    # rasusa downsampling reads to specified coverage level
    call rasusa_task.rasusa as rasusa {
      input:
        read1 = metabuli.metabuli_read1_extract,
        samplename = samplename,
        genome_length = select_first([genome_length])
    }
  }
  # extracted/filtered clean read quality check.
  call nanoplot_task.nanoplot as nanoplot_clean {
    input:
      read1 = select_first([rasusa.read1_subsampled, metabuli.metabuli_read1_extract]),
      samplename = samplename,
      est_genome_length = select_first([genome_length, 12500])
  }
  # check for minimum number of reads, basepairs, coverage, etc
  if (! skip_screen) {
    call screen_task.check_reads_se as clean_check_reads {
      input:
        read1 = select_first([rasusa.read1_subsampled, metabuli.metabuli_read1_extract]),
        workflow_series = "theiaviral",
        expected_genome_length = select_first([genome_length, 12500]), # default to 12500 if genome_length not provided
        skip_mash = true
    }
  }
  if (select_first([clean_check_reads.read_screen, ""]) == "PASS" || skip_screen) {
    # run de novo if no reference genome is provided so we can select a reference
    if (! defined(reference_fasta)) {
      if (call_raven) {
        # de novo assembly with raven
        call raven_task.raven {
          input:
            read1 = select_first([rasusa.read1_subsampled, metabuli.metabuli_read1_extract]),
            samplename = samplename
        }
      }
      if (select_first([raven.raven_status, "FAIL"]) == "FAIL") {
        call flye_task.flye {
          input:
            read1 = select_first([rasusa.read1_subsampled, metabuli.metabuli_read1_extract]),
            samplename = samplename,
            uneven_coverage_mode = true
        }
      }
      # fail gracefully if both assemblies fail
      if (select_first([flye.flye_status, raven.raven_status, "FAIL"]) == "PASS") {
        # quality control metrics for de novo assembly (ie. completeness, viral gene count, contamination)
        call checkv_task.checkv as checkv_denovo {
          input:
            assembly = select_first([flye.assembly_fasta, raven.assembly_fasta]),
            samplename = samplename
        }
        # quality control metrics for de novo assembly (ie. contigs, n50, GC content, genome length)
        call quast_task.quast as quast_denovo {
          input:
            assembly = select_first([flye.assembly_fasta, raven.assembly_fasta]),
            samplename = samplename
        }
      }
    }
    if (defined(reference_fasta) || select_first([flye.flye_status, raven.raven_status, "FAIL"]) == "PASS") {
      # ANI-based reference genome selection
      call skani_task.skani as skani {
        input:
          assembly_fasta = select_first([reference_fasta, flye.assembly_fasta, raven.assembly_fasta]),
          samplename = samplename
      }
      if (defined(reference_fasta) || skani.skani_status == "PASS") {
        # align assembly to reference genome
        call minimap2_task.minimap2 as minimap2 {
          input:
            query1 = select_first([rasusa.read1_subsampled, metabuli.metabuli_read1_extract]),
            reference = select_first([reference_fasta, skani.skani_reference_assembly]),
            samplename = samplename,
            mode = "map-ont",
            output_sam = true,
            long_read_flags = true
        }
        # generate bam file from sam output
        call parse_mapping_task.sam_to_sorted_bam as parse_mapping {
          input:
            sam = minimap2.minimap2_out,
            samplename = samplename,
            min_qual = min_map_quality
        }
        # quality control metrics for reads mapping to reference (ie. coverage, depth, base/map quality)
        call assembly_metrics_task.stats_n_coverage as read_mapping_stats {
          input:
            bamfile = parse_mapping.bam,
            samplename = samplename
        }
        # Index the reference genome for Clair3
        call fasta_utilities_task.samtools_faidx as fasta_utilities{
          input:
            fasta = select_first([reference_fasta, skani.skani_reference_assembly])
        }
        # variant calling with Clair3
        call clair3_task.clair3_variants as clair3 {
          input:
            alignment_bam_file = parse_mapping.bam,
            alignment_bam_file_index = parse_mapping.bai,
            reference_genome_file = select_first([reference_fasta, skani.skani_reference_assembly]),
            reference_genome_file_index = fasta_utilities.fai,
            sequencing_platform = "ont",
            enable_long_indel = true,
            samplename = samplename
        }
        # mask low coverage regions with Ns
        call parse_mapping_task.mask_low_coverage {
          input:
            bam = parse_mapping.bam,
            bai = parse_mapping.bai,
            reference_fasta = select_first([reference_fasta, skani.skani_reference_assembly]),
            min_depth = min_depth
        }
        # create consensus genome based on variant calls
        call bcftools_consensus_task.bcftools_consensus as bcftools_consensus {
          input:
            reference_fasta = mask_low_coverage.mask_reference_fasta,
            input_vcf = clair3.clair3_variants_vcf,
            min_depth = min_depth,
            min_freq = min_allele_freq,
            samplename = samplename
        }
        # quality control metrics for consensus (ie. number of bases, degenerate bases, genome length)
        call consensus_qc_task.consensus_qc as consensus_qc {
          input:
            assembly_fasta = bcftools_consensus.assembly_fasta,
            reference_genome = select_first([reference_fasta, skani.skani_reference_assembly])
        }
        # quality control metrics for consensus (ie. completeness, viral gene count, contamination)
        call checkv_task.checkv as checkv_consensus {
          input:
            assembly = bcftools_consensus.assembly_fasta,
            samplename = samplename
        }
        # run morgana magic for classification
        call morgana_magic_wf.morgana_magic {
          input:
            read1 = select_first([rasusa.read1_subsampled, metabuli.metabuli_read1_extract]),
            samplename = samplename,
            assembly_fasta = select_first([bcftools_consensus.assembly_fasta]),
            taxon_name = ete4_identify.raw_taxon_id,
            seq_method = "nanopore",
            number_ATCG = consensus_qc.number_ATCG,
            workflow_type = "theiaviral",
            reference_gene_locations_bed = reference_gene_locations_bed
        }
      }
    }
  }
  output {
    # versioning outputs
    String theiaviral_ont_version = version_capture.phb_version
    String theiaviral_ont_date = version_capture.date
    # ncbi datasets - taxon identification
    String ncbi_taxon_id = ete4_identify.taxon_id
    String ncbi_taxon_name = ete4_identify.taxon_name
    String ncbi_read_extraction_rank = ete4_identify.taxon_rank
    String ete4_version = ete4_identify.ete4_version
    String ete4_docker = ete4_identify.ete4_docker
    # host decontamination outputs
    File? dehost_wf_dehost_read1 = host_decontaminate.dehost_read1
    String? dehost_wf_host_accession = host_decontaminate.host_genome_accession
    File? dehost_wf_host_fasta = host_decontaminate.host_genome_fasta
    File? dehost_wf_host_mapped_bam = host_decontaminate.host_mapped_sorted_bam
    File? dehost_wf_host_mapped_bai = host_decontaminate.host_mapped_sorted_bai
    File? dehost_wf_host_mapping_stats = host_decontaminate.host_mapping_stats
    File? dehost_wf_host_mapping_cov_hist = host_decontaminate.host_mapping_cov_hist
    File? dehost_wf_host_flagstat = host_decontaminate.host_flagstat
    Float? dehost_wf_host_mapping_coverage = host_decontaminate.host_mapping_coverage
    Float? dehost_wf_host_mapping_mean_depth = host_decontaminate.host_mapping_mean_depth
    Float? dehost_wf_host_percent_mapped_reads = host_decontaminate.host_percent_mapped_reads
    File? dehost_wf_host_mapping_metrics = host_decontaminate.host_mapping_metrics
    # raw read quality control
    File nanoplot_html_raw = nanoplot_raw.nanoplot_html
    File nanoplot_tsv_raw = nanoplot_raw.nanoplot_tsv
    Int nanoplot_num_reads_raw1 = nanoplot_raw.num_reads
    Float nanoplot_r1_median_readlength_raw = nanoplot_raw.median_readlength
    Float nanoplot_r1_mean_readlength_raw = nanoplot_raw.mean_readlength
    Float nanoplot_r1_stdev_readlength_raw = nanoplot_raw.stdev_readlength
    Float nanoplot_r1_n50_raw = nanoplot_raw.n50
    Float nanoplot_r1_mean_q_raw = nanoplot_raw.mean_q
    Float nanoplot_r1_median_q_raw = nanoplot_raw.median_q
    # porechop outputs - adapter trimming
    File? porechop_trimmed_read1 = porechop.trimmed_reads
    String? porechop_version = porechop.porechop_version
    # nanoq outputs - read filtering
    File nanoq_filtered_read1 = nanoq.filtered_read1
    String nanoq_version = nanoq.version
    # scrubbed reads
    File ncbi_scrub_read1_dehosted = ncbi_scrub_se.read1_dehosted
    Int ncbi_scrub_human_spots_removed = ncbi_scrub_se.human_spots_removed
    String ncbi_scrub_docker = ncbi_scrub_se.ncbi_scrub_docker
    # metabuli outputs - taxonomic classification and read extraction
    File? metabuli_report = metabuli.metabuli_report
    File? metabuli_classified = metabuli.metabuli_classified
    File? metabuli_read1_extract = metabuli.metabuli_read1_extract
    File? metabuli_krona_report = metabuli.metabuli_krona_report
    String? metabuli_database = metabuli.metabuli_database
    String? metabuli_version = metabuli.metabuli_version
    String? metabuli_docker = metabuli.metabuli_docker
    # rasusa outputs - downsampled reads
    File? rasusa_read1_subsampled = rasusa.read1_subsampled
    File? rasusa_read2_subsampled = rasusa.read2_subsampled
    String? rasusa_version = rasusa.rasusa_version
    # clean read quality control
    File? nanoplot_html_clean = nanoplot_clean.nanoplot_html
    File? nanoplot_tsv_clean = nanoplot_clean.nanoplot_tsv
    Int? nanoplot_num_reads_clean1 = nanoplot_clean.num_reads
    Float? nanoplot_r1_median_readlength_clean = nanoplot_clean.median_readlength
    Float? nanoplot_r1_mean_readlength_clean = nanoplot_clean.mean_readlength
    Float? nanoplot_r1_stdev_readlength_clean = nanoplot_clean.stdev_readlength
    Float? nanoplot_r1_n50_clean = nanoplot_clean.n50
    Float? nanoplot_r1_mean_q_clean = nanoplot_clean.mean_q
    Float? nanoplot_r1_median_q_clean = nanoplot_clean.median_q
    # clean read screen outputs
    String? read_screen_clean = clean_check_reads.read_screen
    File? read_screen_clean_tsv = clean_check_reads.read_screen_tsv
    # raven outputs - denovo genome assembly
    File? assembly_denovo_fasta = select_first([flye.assembly_fasta, raven.assembly_fasta, "gs://theiagen-public-resources-rp/empty_files/empty.fasta"])
    String? raven_denovo_status = raven.raven_status
    String? raven_denovo_version = raven.raven_version
    String? raven_denovo_docker = raven.raven_docker
    String? flye_denovo_status = flye.flye_status
    String? flye_denovo_version = flye.flye_version
    String? flye_denovo_docker = flye.flye_docker
    File? flye_denovo_info = flye.assembly_info
    # checkv_denovo outputs - denovo assembly quality control
    File? checkv_denovo_summary = checkv_denovo.checkv_summary
    File? checkv_denovo_contamination = checkv_denovo.checkv_contamination
    Float? checkv_denovo_weighted_contamination = checkv_denovo.weighted_contamination
    Float? checkv_denovo_weighted_completeness = checkv_denovo.weighted_completeness
    Int? checkv_denovo_total_genes = checkv_denovo.total_genes
    String? checkv_denovo_version = checkv_denovo.checkv_version
    # quast_denovo outputs - denovo assembly quality control
    File? quast_denovo_report = quast_denovo.quast_report
    Int? quast_denovo_genome_length = quast_denovo.genome_length
    Int? quast_denovo_number_contigs = quast_denovo.number_contigs
    Int? quast_denovo_n50_value = quast_denovo.n50_value
    Int? quast_denovo_largest_contig = quast_denovo.largest_contig
    Float? quast_denovo_gc_percent = quast_denovo.gc_percent
    Float? quast_denovo_uncalled_bases = quast_denovo.uncalled_bases
    String? quast_denovo_version = quast_denovo.version
    String? quast_denovo_docker = quast_denovo.quast_docker
    # skani outputs - ANI-based reference genome selection
    File? skani_report = skani.skani_report
    String? skani_top_accession = skani.skani_top_accession
    String? skani_reference_taxon = skani.skani_reference_taxon
    Float? skani_top_score = skani.skani_top_score
    Float? skani_top_ani = skani.skani_top_ani
    Float? skani_top_query_coverage = skani.skani_top_query_coverage
    String? skani_warning = skani.skani_warning
    String? skani_status = skani.skani_status
    String? skani_database = skani.skani_database
    String? skani_version = skani.skani_version
    String? skani_docker = skani.skani_docker
    File? skani_reference_assembly = skani.skani_reference_assembly
    # minimap2 outputs - reads aligned to best reference
    File? minimap2_out = minimap2.minimap2_out
    String? minimap2_version = minimap2.minimap2_version
    String? minimap2_docker = minimap2.minimap2_docker
    # parse_mapping outputs - sam to sorted bam conversion
    File? assembly_to_ref_bam = parse_mapping.bam
    File? assembly_to_ref_bai = parse_mapping.bai
    String? parse_mapping_samtools_version = parse_mapping.samtools_version
    String? parse_mapping_samtools_docker = parse_mapping.samtools_docker
    # assembly_metrics outputs - read mapping quality control
    File? read_mapping_report = read_mapping_stats.metrics_txt
    File? read_mapping_statistics = read_mapping_stats.stats
    File? read_mapping_cov_hist = read_mapping_stats.cov_hist
    File? read_mapping_cov_stats = read_mapping_stats.cov_stats
    File? read_mapping_flagstat = read_mapping_stats.flagstat
    Float? read_mapping_coverage = read_mapping_stats.coverage
    Float? read_mapping_depth = read_mapping_stats.depth
    Float? read_mapping_meanbaseq = read_mapping_stats.meanbaseq
    Float? read_mapping_meanmapq = read_mapping_stats.meanmapq
    Float? read_mapping_percentage_mapped_reads = read_mapping_stats.percentage_mapped_reads
    String? read_mapping_date = read_mapping_stats.date
    String? read_mapping_samtools_version = read_mapping_stats.samtools_version
    # fasta_utilities outputs - samtools faidx reference genome
    File? fasta_utilities_fai = fasta_utilities.fai
    String? fasta_utilities_samtools_version = fasta_utilities.samtools_version
    String? fasta_utilities_samtools_docker = fasta_utilities.samtools_docker
    # clair3 outputs - variant calling
    File? clair3_vcf = clair3.clair3_variants_vcf
    File? clair3_gvcf = clair3.clair3_variants_gvcf
    String? clair3_model = clair3.clair3_model_used
    String? clair3_version = clair3.clair3_version
    String? clair3_docker = clair3.clair3_variants_docker_image
    # coverage_mask outputs - low coverage regions
    File? mask_low_coverage_bed = mask_low_coverage.low_coverage_regions_bed
    File? mask_low_coverage_all_coverage_bed = mask_low_coverage.all_coverage_regions_bed
    File? mask_low_coverage_reference_fasta = mask_low_coverage.mask_reference_fasta
    String? mask_low_coverage_bedtools_version = mask_low_coverage.bedtools_version
    String? mask_low_coverage_bedtools_docker = mask_low_coverage.bedtools_docker
    # bcftools_consensus outputs - consensus genome
    File? assembly_consensus_fasta = bcftools_consensus.assembly_fasta
    File? bcftools_filtered_vcf = bcftools_consensus.bcftools_filtered_vcf
    String? bcftools_version = bcftools_consensus.bcftools_version
    String? bcftools_docker = bcftools_consensus.bcftools_docker
    # consensus assembly statistics
    Int? consensus_qc_number_N = consensus_qc.number_N
    Int? consensus_qc_assembly_length_unambiguous = consensus_qc.number_ATCG
    Int? consensus_qc_number_Degenerate = consensus_qc.number_Degenerate
    Int? consensus_qc_number_Total = consensus_qc.number_Total
    Float? consensus_qc_percent_reference_coverage = consensus_qc.percent_reference_coverage
    # checkv_consensus outputs - consensus assembly quality control
    File? checkv_consensus_summary = checkv_consensus.checkv_summary
    File? checkv_consensus_contamination = checkv_consensus.checkv_contamination
    Float? checkv_consensus_weighted_contamination = checkv_consensus.weighted_contamination
    Float? checkv_consensus_weighted_completeness = checkv_consensus.weighted_completeness
    Int? checkv_consensus_total_genes = checkv_consensus.total_genes
    String? checkv_consensus_version = checkv_consensus.checkv_version
    # morgana magic outputs
    String? morgana_magic_organism = morgana_magic.organism
    # VADR outputs
    File? vadr_alerts_list = morgana_magic.vadr_alerts_list
    String? vadr_num_alerts = morgana_magic.vadr_num_alerts
    File? vadr_feature_tbl_pass = morgana_magic.vadr_feature_tbl_pass
    File? vadr_feature_tbl_fail = morgana_magic.vadr_feature_tbl_fail
    File? vadr_classification_summary_file = morgana_magic.vadr_classification_summary_file
    File? vadr_all_outputs_tar_gz = morgana_magic.vadr_all_outputs_tar_gz
    String? vadr_docker = morgana_magic.vadr_docker
    File? vadr_fastas_zip_archive = morgana_magic.vadr_fastas_zip_archive
    # Pangolin outputs
    String? pango_lineage = morgana_magic.pango_lineage
    String? pango_lineage_expanded = morgana_magic.pango_lineage_expanded
    String? pangolin_conflicts = morgana_magic.pangolin_conflicts
    String? pangolin_notes = morgana_magic.pangolin_notes
    String? pangolin_assignment_version = morgana_magic.pangolin_assignment_version
    File? pango_lineage_report = morgana_magic.pango_lineage_report
    String? pangolin_docker = morgana_magic.pangolin_docker
    String? pangolin_versions = morgana_magic.pangolin_versions
    # Nextclade outputs for all organisms
    String? nextclade_version = morgana_magic.nextclade_version
    String? nextclade_docker = morgana_magic.nextclade_docker
    String? nextclade_ds_tag = morgana_magic.nextclade_ds_tag
    File? nextclade_json = morgana_magic.nextclade_json
    File? auspice_json = morgana_magic.auspice_json
    File? nextclade_tsv = morgana_magic.nextclade_tsv
    String? nextclade_aa_subs = morgana_magic.nextclade_aa_subs
    String? nextclade_aa_dels = morgana_magic.nextclade_aa_dels
    String? nextclade_clade = morgana_magic.nextclade_clade
    String? nextclade_lineage = morgana_magic.nextclade_lineage
    String? nextclade_qc = morgana_magic.nextclade_qc
    # Nextclade outputs for Rabies
    File? nextclade_json_rabies = morgana_magic.nextclade_json_rabies
    File? auspice_json_rabies = morgana_magic.auspice_json_rabies
    File? nextclade_tsv_rabies = morgana_magic.nextclade_tsv_rabies
    String? nextclade_aa_subs_rabies = morgana_magic.nextclade_aa_subs_rabies
    String? nextclade_aa_dels_rabies = morgana_magic.nextclade_aa_dels_rabies
    String? nextclade_clade_rabies = morgana_magic.nextclade_clade_rabies
    String? nextclade_lineage_rabies = morgana_magic.nextclade_lineage_rabies
    String? nextclade_qc_rabies = morgana_magic.nextclade_qc_rabies
    # Nextclade outputs for flu HA
    File? nextclade_json_flu_ha = morgana_magic.nextclade_json_flu_ha
    File? auspice_json_flu_ha = morgana_magic.auspice_json_flu_ha
    File? nextclade_tsv_flu_ha = morgana_magic.nextclade_tsv_flu_ha
    String? nextclade_ds_tag_flu_ha = morgana_magic.nextclade_ds_tag_flu_ha
    String? nextclade_aa_subs_flu_ha = morgana_magic.nextclade_aa_subs_flu_ha
    String? nextclade_aa_dels_flu_ha = morgana_magic.nextclade_aa_dels_flu_ha
    String? nextclade_clade_flu_ha = morgana_magic.nextclade_clade_flu_ha
    String? nextclade_qc_flu_ha = morgana_magic.nextclade_qc_flu_ha
    # Nextclade outputs for flu NA
    File? nextclade_json_flu_na = morgana_magic.nextclade_json_flu_na
    File? auspice_json_flu_na = morgana_magic.auspice_json_flu_na
    File? nextclade_tsv_flu_na = morgana_magic.nextclade_tsv_flu_na
    String? nextclade_ds_tag_flu_na = morgana_magic.nextclade_ds_tag_flu_na
    String? nextclade_aa_subs_flu_na = morgana_magic.nextclade_aa_subs_flu_na
    String? nextclade_aa_dels_flu_na = morgana_magic.nextclade_aa_dels_flu_na
    String? nextclade_clade_flu_na = morgana_magic.nextclade_clade_flu_na
    String? nextclade_qc_flu_na = morgana_magic.nextclade_qc_flu_na
    # Flu IRMA Outputs
    String? irma_version = morgana_magic.irma_version
    String? irma_docker = morgana_magic.irma_docker
    String? irma_type = morgana_magic.irma_type
    String? irma_subtype = morgana_magic.irma_subtype
    String? irma_subtype_notes = morgana_magic.irma_subtype_notes
    # Flu GenoFLU Outputs
    String? genoflu_version = morgana_magic.genoflu_version
    String? genoflu_genotype = morgana_magic.genoflu_genotype
    String? genoflu_all_segments = morgana_magic.genoflu_all_segments
    File? genoflu_output_tsv = morgana_magic.genoflu_output_tsv
    # Flu Abricate Outputs
    String? abricate_flu_type = morgana_magic.abricate_flu_type
    String? abricate_flu_subtype =  morgana_magic.abricate_flu_subtype
    File? abricate_flu_results = morgana_magic.abricate_flu_results
    String? abricate_flu_database =  morgana_magic.abricate_flu_database
    String? abricate_flu_version = morgana_magic.abricate_flu_version
    # HIV Quasitools Outputs
    String? quasitools_version = morgana_magic.quasitools_version
    String? quasitools_date = morgana_magic.quasitools_date
    File? quasitools_coverage_file = morgana_magic.quasitools_coverage_file
    File? quasitools_dr_report = morgana_magic.quasitools_dr_report
    File? quasitools_hydra_vcf = morgana_magic.quasitools_hydra_vcf
    File? quasitools_mutations_report = morgana_magic.quasitools_mutations_report
  }
}