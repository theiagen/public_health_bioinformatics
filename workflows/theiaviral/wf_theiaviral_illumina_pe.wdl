version 1.0

import "../../tasks/task_versioning.wdl" as versioning_task
import "../../tasks/quality_control/comparisons/task_screen.wdl" as read_screen_task
import "../../tasks/quality_control/basic_statistics/task_fastq_scan.wdl" as fastq_scan_task
import "../../tasks/quality_control/read_filtering/task_bbduk.wdl" as bbduk_task
import "../../tasks/quality_control/read_filtering/task_fastp.wdl" as fastp_task
import "../../tasks/quality_control/read_filtering/task_ncbi_scrub.wdl" as ncbi_scrub_task
import "../../tasks/taxon_id/contamination/task_kraken2.wdl" as kraken2_task
import "../../tasks/utilities/file_handling/task_cat_lanes.wdl" as cat_lanes_task
import "../../tasks/taxon_id/task_krakentools.wdl" as krakentools_task
import "../../tasks/utilities/task_rasusa.wdl" as rasusa_task
import "../../tasks/assembly/task_spades.wdl" as spades_task
import "../../tasks/assembly/task_megahit.wdl" as megahit_task
import "../../tasks/quality_control/advanced_metrics/task_checkv.wdl" as checkv_task
import "../../tasks/quality_control/basic_statistics/task_quast.wdl" as quast_task
import "../../tasks/taxon_id/task_skani.wdl" as skani_task
import "../../tasks/taxon_id/task_ete4_taxon_id.wdl" as identify_taxon_id_task
import "../../tasks/utilities/task_datasets_genome_length.wdl" as genome_length_task
import "../../tasks/alignment/task_bwa.wdl" as bwa_task
import "../../tasks/assembly/task_ivar_consensus.wdl" as ivar_consensus_task
import "../../tasks/gene_typing/variant_detection/task_ivar_variant_call.wdl" as variant_call_task
import "../../tasks/quality_control/basic_statistics/task_assembly_metrics.wdl" as assembly_metrics_task
import "../../tasks/quality_control/basic_statistics/task_consensus_qc.wdl" as consensus_qc_task
import "../utilities/wf_host_decontaminate.wdl" as host_decontaminate_wf
import "../utilities/wf_morgana_magic.wdl" as morgana_magic_wf

workflow theiaviral_illumina_pe {
  meta {
    description: "De novo assembly, dynamic reference selection, and reference-based consensus calling for viral metagenomic/genomic data generated on Illumina paired-end NGS platforms."
  }
  input {
    File read1
    File read2
    String taxon # taxon id OR organism name
    String read_extraction_rank = "family"
    String samplename
    String? host # host to dehost reads, if provided
    File kraken_db = "gs://theiagen-public-resources-rp/reference_data/databases/kraken2/kraken2_humanGRCh38_viralRefSeq_20240828.tar.gz"
    File? checkv_db 
    File? skani_db
    Boolean skip_screen = false # if false, run clean read screening
    Boolean skip_qc = false # If false, run read quality control
    Boolean call_metaviralspades = true # if false, move to megahit immediately
    Boolean skip_rasusa = true 
    File? reference_fasta # optional, if provided, will be used instead of dynamic reference selection
    File? reference_gene_locations_bed # optional, if provided will be used for coverage calculations
    Boolean extract_unclassified = true # if true, unclassified reads will be extracted from kraken2 output
    Int min_depth = 10
    Int min_map_quality = 20
    Float min_allele_freq = 0.6
    # rasusa downsampling inputs
    Int? genome_length
  }
  # get the PHB version
  call versioning_task.version_capture {
    input:
  }
  if (!skip_qc) {
    # get the taxon id
    call identify_taxon_id_task.ete4_taxon_id as ete4_identify {
      input:
        taxon = taxon,
        rank = read_extraction_rank
    }
    if (! defined(genome_length)) {
      # get average genome length for the taxon
      call genome_length_task.datasets_genome_length as est_genome_length {
        input:
          taxon = select_first([ete4_identify.raw_taxon_id, taxon]),
          use_ncbi_virus = true,
          complete = true,
          refseq = true
      }
    }
    # read QC, classification, extraction, and trimming
    call fastq_scan_task.fastq_scan_pe as fastq_scan_raw {
      input:
        read1 = read1,
        read2 = read2
    }
    # human read scrubbing
    call ncbi_scrub_task.ncbi_scrub_pe {
      input:
        samplename = samplename,
        read1 = read1,
        read2 = read2
    }
    # adapter + read trimming
    call fastp_task.fastp_pe as fastp {
      input:
        samplename = samplename,
        read1 = ncbi_scrub_pe.read1_dehosted,
        read2 = ncbi_scrub_pe.read2_dehosted
    }
    call bbduk_task.bbduk {
      input:
        samplename = samplename,
        read1_trimmed = fastp.read1_trimmed,
        read2_trimmed = fastp.read2_trimmed
    }
    # host read decontamination
    if (defined(host)) {
      call host_decontaminate_wf.host_decontaminate {
        input:
          samplename = samplename,
          read1 = bbduk.read1_clean,
          read2 = bbduk.read2_clean,
          host = select_first([host])
      }
    }
    # taxon-based read extraction
    call kraken2_task.kraken2_standalone {
      input:
        samplename = samplename,
        read1 = select_first([host_decontaminate.dehost_read1, bbduk.read1_clean]),
        read2 = select_first([host_decontaminate.dehost_read2, bbduk.read2_clean]),
        kraken2_db = kraken_db
    }
    call krakentools_task.extract_kraken_reads as kraken2_extract {
      input:
        read1 = kraken2_standalone.kraken2_classified_read1,
        read2 = select_first([kraken2_standalone.kraken2_classified_read2]),
        taxon_id = ete4_identify.taxon_id,
        kraken2_output = kraken2_standalone.kraken2_classified_report,
        kraken2_report = kraken2_standalone.kraken2_report
    }
    # inclusion of unclassified reads
    if (extract_unclassified) {
      call cat_lanes_task.cat_lanes {
        input:
          samplename = samplename,
          read1_lane1 = kraken2_standalone.kraken2_unclassified_read1,
          read1_lane2 = select_first([kraken2_extract.extracted_read1]),
          read2_lane1 = kraken2_standalone.kraken2_unclassified_read2,
          read2_lane2 = select_first([kraken2_extract.extracted_read2])
      }
    }
    if (! skip_rasusa) {
      # downsample reads to a specific coverage
      call rasusa_task.rasusa as rasusa {
        input:
          read1 = select_first([cat_lanes.read1_concatenated, kraken2_extract.extracted_read1, read1]),
          read2 = select_first([cat_lanes.read2_concatenated, kraken2_extract.extracted_read2, read2]),
          samplename = samplename,
          genome_length = select_first([est_genome_length.avg_genome_length, genome_length])
      }
    }
    # clean read stat gathering
    call fastq_scan_task.fastq_scan_pe as fastq_scan_clean {
      input:
        read1 = select_first([rasusa.read1_subsampled, cat_lanes.read1_concatenated, kraken2_extract.extracted_read1]),
        read2 = select_first([rasusa.read2_subsampled, cat_lanes.read2_concatenated, kraken2_extract.extracted_read2])
    }
  }
  # clean read screening
  if (! skip_screen) {
    call read_screen_task.check_reads as clean_check_reads {
      input:
        read1 = select_first([rasusa.read1_subsampled, cat_lanes.read1_concatenated, kraken2_extract.extracted_read1, read1]),
        read2 = select_first([rasusa.read2_subsampled, cat_lanes.read2_concatenated, kraken2_extract.extracted_read2, read2]),
        workflow_series = "theiaviral",
        expected_genome_length = select_first([est_genome_length.avg_genome_length, genome_length])
    }
  }
  if (select_first([clean_check_reads.read_screen, ""]) == "PASS" || skip_screen) {
    # run de novo if no reference genome is provided so we can select a reference
    if (! defined(reference_fasta)) {
      # de novo assembly - prioritize metaviralspades
      if (call_metaviralspades) {
        call spades_task.spades {
          input:
            read1 = select_first([rasusa.read1_subsampled, cat_lanes.read1_concatenated, kraken2_extract.extracted_read1, read1]),
            read2 = select_first([rasusa.read2_subsampled, cat_lanes.read2_concatenated, kraken2_extract.extracted_read2, read2]),
            samplename = samplename,
            spades_type = "metaviral"
        }
      }
      # fallback to megahit if metaviralspades fails to identify a complete virus
      if (select_first([spades.spades_status, "FAIL"]) == "FAIL") {
        call megahit_task.megahit {
          input:
            read1 = select_first([rasusa.read1_subsampled, cat_lanes.read1_concatenated, kraken2_extract.extracted_read1, read1]),
            read2 = select_first([rasusa.read2_subsampled, cat_lanes.read2_concatenated, kraken2_extract.extracted_read2, read2]),
            samplename = samplename
        }
      }
      # fail gracefully if both assemblies fail
      if (select_first([megahit.megahit_status, spades.spades_status, "FAIL"]) == "PASS") {
        # quality control metrics for de novo assembly (ie. completeness, viral gene count, contamination)
        call checkv_task.checkv as checkv_denovo {
          input:
            assembly = select_first([spades.assembly_fasta, megahit.assembly_fasta]),
            samplename = samplename,
            checkv_db = checkv_db
        }
        # quality control metrics for de novo assembly (ie. contigs, n50, GC content, genome length)
        call quast_task.quast as quast_denovo {
          input:
            assembly = select_first([spades.assembly_fasta, megahit.assembly_fasta]),
            samplename = samplename,
            min_contig_length = 0
        }
      }
    }
    if (defined(reference_fasta) || select_first([megahit.megahit_status, spades.spades_status, "FAIL"]) == "PASS") {
      # ANI-based reference genome selection
      call skani_task.skani as skani {
        input:
          assembly_fasta = select_first([reference_fasta, spades.assembly_fasta, megahit.assembly_fasta]),
          samplename = samplename,
          skani_db = skani_db
      }
      if (defined(reference_fasta) || skani.skani_status == "PASS") {
        # align reads to reference
        call bwa_task.bwa {
          input:
            samplename = samplename,
            read1 = select_first([rasusa.read1_subsampled, cat_lanes.read1_concatenated, kraken2_extract.extracted_read1, read1]),
            read2 = select_first([rasusa.read2_subsampled, cat_lanes.read2_concatenated, kraken2_extract.extracted_read2, read2]),
            reference_genome = select_first([reference_fasta, skani.skani_reference_assembly])
        }
        # consensus calling via ivar
        call ivar_consensus_task.consensus {
          input:
            bamfile = bwa.sorted_bam,
            samplename = samplename,
            reference_genome = select_first([reference_fasta, skani.skani_reference_assembly]),
            min_qual = min_map_quality,
            consensus_min_depth = min_depth,
            consensus_min_freq = min_allele_freq,
            all_positions = true
        }
        # variant calling via ivar
        call variant_call_task.variant_call as ivar_variants {
          input:
            mpileup = consensus.sample_mpileup,
            samplename = samplename,
            reference_genome = select_first([reference_fasta, skani.skani_reference_assembly]),
            min_qual = min_map_quality,
            organism = "",
            variant_min_freq = min_allele_freq,
            variant_min_depth = min_depth
        }
        # quality control metrics for reads mapping to reference (ie. coverage, depth, base/map quality)
        call assembly_metrics_task.stats_n_coverage as read_mapping_stats {
          input:
            bamfile = bwa.sorted_bam,
            samplename = samplename,
            read1 = select_first([rasusa.read1_subsampled, cat_lanes.read1_concatenated, kraken2_extract.extracted_read1, read1]),
            read2 = select_first([rasusa.read2_subsampled, cat_lanes.read2_concatenated, kraken2_extract.extracted_read2, read2]),
        }
        # quality control metrics for consensus (ie. number of bases, degenerate bases, genome length)
        call consensus_qc_task.consensus_qc as consensus_qc {
          input:
            assembly_fasta = consensus.consensus_seq,
            reference_genome = select_first([reference_fasta, skani.skani_reference_assembly]),
            genome_length = select_first([est_genome_length.avg_genome_length, genome_length])
        }
        # quality control metrics for consensus (ie. completeness, viral gene count, contamination)
        call checkv_task.checkv as checkv_consensus {
          input:
            assembly = consensus.consensus_seq,
            samplename = samplename,
            checkv_db = checkv_db
        }
        # run morgana magic for classification
        call morgana_magic_wf.morgana_magic {
          input:
            samplename = samplename,
            read1 = select_first([rasusa.read1_subsampled, cat_lanes.read1_concatenated, kraken2_extract.extracted_read1, read1]),
            read2 = select_first([rasusa.read2_subsampled, cat_lanes.read2_concatenated, kraken2_extract.extracted_read2, read2]),
            assembly_fasta = select_first([consensus.consensus_seq]),
            taxon_name = select_first([ete4_identify.raw_taxon_id, taxon]),
            seq_method = "illumina_pe",
            number_ATCG = consensus_qc.number_ATCG,
            workflow_type = "theiaviral",
            reference_gene_locations_bed = reference_gene_locations_bed
        }
      }
    }
  }
  output {
    # versioning outputs
    String theiaviral_illumina_pe_version = version_capture.phb_version
    String theiaviral_illumina_pe_date = version_capture.date
    # ete4 - taxon identification
    String? ncbi_taxon_id = ete4_identify.taxon_id
    String? ncbi_taxon_name = ete4_identify.taxon_name
    String? ncbi_read_extraction_rank = ete4_identify.taxon_rank
    String? ete4_version = ete4_identify.ete4_version
    String? ete4_docker = ete4_identify.ete4_docker
    # NCBI datasets genome length estimation
    String? taxon_avg_genome_length = est_genome_length.avg_genome_length
    String? datasets_genome_length_docker = est_genome_length.ncbi_datasets_docker
    String? datasets_genome_length_version = est_genome_length.ncbi_datasets_version
    # raw read quality control
    Int? fastq_scan_num_reads_raw1 = fastq_scan_raw.read1_seq
    Int? fastq_scan_num_reads_raw2 = fastq_scan_raw.read2_seq
    String? fastq_scan_raw_pairs = fastq_scan_raw.read_pairs
    String? fastq_scan_version = fastq_scan_raw.version
    String? fastq_scan_docker = fastq_scan_raw.fastq_scan_docker
    File? fastq_scan_raw1_json = fastq_scan_raw.read1_fastq_scan_json
    File? fastq_scan_raw2_json = fastq_scan_raw.read2_fastq_scan_json
    # NCBI scrubbing
    File? ncbi_scrub_read1_dehosted = ncbi_scrub_pe.read1_dehosted
    File? ncbi_scrub_read2_dehosted = ncbi_scrub_pe.read2_dehosted
    Int? ncbi_scrub_human_spots_removed = ncbi_scrub_pe.human_spots_removed
    String? ncbi_scrub_docker = ncbi_scrub_pe.ncbi_scrub_docker
    # host decontamination outputs
    File? dehost_wf_dehost_read1 = host_decontaminate.dehost_read1
    File? dehost_wf_dehost_read2 = host_decontaminate.dehost_read2
    String? dehost_wf_host_accession = host_decontaminate.host_genome_accession
    File? dehost_wf_host_mapped_bam = host_decontaminate.host_mapped_sorted_bam
    File? dehost_wf_host_mapped_bai = host_decontaminate.host_mapped_sorted_bai
    File? dehost_wf_host_fasta = host_decontaminate.host_genome_fasta
    File? dehost_wf_host_mapping_stats = host_decontaminate.host_mapping_stats
    File? dehost_wf_host_mapping_cov_hist = host_decontaminate.host_mapping_cov_hist
    File? dehost_wf_host_flagstat = host_decontaminate.host_flagstat
    Float? dehost_wf_host_mapping_coverage = host_decontaminate.host_mapping_coverage
    Float? dehost_wf_host_mapping_mean_depth = host_decontaminate.host_mapping_mean_depth
    Float? dehost_wf_host_percent_mapped_reads = host_decontaminate.host_percent_mapped_reads
    File? dehost_wf_host_mapping_metrics = host_decontaminate.host_mapping_metrics
    # trimming outputs - adapter trimming
    String? fastp_version = fastp.version
    File? fastp_html_report = fastp.fastp_stats
    # bbduk outputs
    String? bbduk_docker = bbduk.bbduk_docker
    File? bbduk_read1_clean = bbduk.read1_clean
    File? bbduk_read2_clean = bbduk.read2_clean
    # kraken2 outputs - taxonomic classification and read extraction
    String? kraken2_version = kraken2_standalone.kraken2_version
    String? kraken2_docker = kraken2_standalone.kraken2_docker
    String? kraken2_database = kraken2_standalone.kraken2_database
    String? kraken2_report = kraken2_standalone.kraken2_report
    File? kraken2_extracted_read1 = kraken2_extract.extracted_read1
    File? kraken2_extracted_read2 = kraken2_extract.extracted_read2
    String? kraken2_extraction_status = kraken2_extract.status
    # clean read quality control
    Int? fastq_scan_num_reads_clean1 = fastq_scan_clean.read1_seq
    Int? fastq_scan_num_reads_clean2 = fastq_scan_clean.read2_seq
    String? fastq_scan_clean_pairs = fastq_scan_clean.read_pairs
    File? fastq_scan_clean1_json = fastq_scan_clean.read1_fastq_scan_json
    File? fastq_scan_clean2_json = fastq_scan_clean.read2_fastq_scan_json
    # subsampled reads
    File? downsampled_read1 = rasusa.read1_subsampled
    File? downsampled_read2 = rasusa.read2_subsampled
    # clean read screening
    String? read_screen_clean = clean_check_reads.read_screen
    File? read_screen_clean_tsv = clean_check_reads.read_screen_tsv
    # denovo genome assembly
    File assembly_denovo_fasta = select_first([spades.assembly_fasta, megahit.assembly_fasta, "gs://theiagen-public-resources-rp/empty_files/empty.fasta"])
    String? metaviralspades_status = spades.spades_status
    String? metaviralspades_version = spades.spades_version
    String? metaviralspades_docker = spades.spades_docker
    String? megahit_status = megahit.megahit_status
    String? megahit_version = megahit.megahit_version
    String? megahit_docker = megahit.megahit_docker
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
    # bwa outputs - reads aligned to best reference
    String? bwa_version = bwa.bwa_version
    String? bwa_samtools_version = bwa.sam_version
    File? bwa_read1_aligned = bwa.read1_aligned
    File? bwa_read2_aligned = bwa.read2_aligned
    File? bwa_sorted_bam = bwa.sorted_bam
    File? bwa_sorted_bai = bwa.sorted_bai
    File? bwa_read1_unaligned = bwa.read1_unaligned
    File? bwa_read2_unaligned = bwa.read2_unaligned
    File? bwa_sorted_bam_unaligned = bwa.sorted_bam_unaligned
    File? bwa_sorted_bam_unaligned_bai = bwa.sorted_bam_unaligned_bai
    # Read mapping stats
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
    # Read Alignment - variant call outputs
    File? ivar_tsv = ivar_variants.sample_variants_tsv
    File? ivar_vcf = ivar_variants.sample_variants_vcf
    String? ivar_variant_proportion_intermediate = ivar_variants.variant_proportion_intermediate
    String? ivar_variant_version = ivar_variants.ivar_version
    # Consensus assembly outputs
    File? assembly_consensus_fasta = consensus.consensus_seq
    Int consensus_n_variant_min_depth = min_depth
    String? ivar_version_consensus = consensus.ivar_version
    # Consensus QC
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