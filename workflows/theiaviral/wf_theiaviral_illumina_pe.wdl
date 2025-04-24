version 1.0

import "../utilities/wf_read_QC_trim_pe.wdl" as read_qc
import "../../tasks/quality_control/comparisons/task_screen.wdl" as read_screen_task
import "../../tasks/utilities/task_rasusa.wdl" as rasusa_task
import "../../tasks/assembly/task_spades.wdl" as spades_task
import "../../tasks/assembly/task_megahit.wdl" as megahit_task
import "../../tasks/quality_control/advanced_metrics/task_checkv.wdl" as checkv_task
import "../../tasks/quality_control/basic_statistics/task_quast.wdl" as quast_task
import "../../tasks/taxon_id/task_skani.wdl" as skani_task
import "../../tasks/utilities/data_import/task_ncbi_datasets.wdl" as ncbi_datasets_task
import "../../tasks/taxon_id/task_identify_taxon_id.wdl" as identify_taxon_id_task
import "../../tasks/alignment/task_bwa.wdl" as bwa_task
import "../../tasks/assembly/task_ivar_consensus.wdl" as ivar_consensus
import "../../tasks/gene_typing/variant_detection/task_ivar_variant_call.wdl" as variant_call_task
import "../../tasks/quality_control/basic_statistics/task_assembly_metrics.wdl" as assembly_metrics_task
import "../../tasks/utilities/data_handling/task_fasta_utilities.wdl" as fasta_utilities_task
import "../../tasks/quality_control/basic_statistics/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/task_versioning.wdl" as versioning


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
    File kraken_db = "gs://theiagen-large-public-files-rp/terra/databases/kraken2/kraken2_humanGRCh38_viralRefSeq_20240828.tar.gz"
    Boolean skip_screen = true # if false, run clean read screening
    Boolean skip_metaviralspades = false # if true, move to megahit immediately
    Boolean skip_rasusa = false
    File? reference_fasta # optional, if provided, will be used instead of dynamic reference selection
    Boolean extract_unclassified = false # if true, unclassified reads will be extracted from kraken2 output
    Int min_depth = 10
    Int min_map_quality = 20
    Float min_allele_freq = 0.6
    # rasusa downsampling inputs
    Int? genome_length
  }
  # get the PHB version
  call versioning.version_capture {
    input:
  }
  # get the taxon id
  call identify_taxon_id_task.identify_taxon_id as ncbi_identify {
    input:
      taxon = taxon,
      rank = read_extraction_rank
  }
  # read QC, classification, extraction, and trimming
  call read_qc.read_QC_trim_pe as read_QC_trim {
    input:
      read1 = read1,
      read2 = read2,
      samplename = samplename,
      taxon_id = ncbi_identify.taxon_id,
      extract_unclassified = extract_unclassified,
      kraken_db = kraken_db,
      workflow_series = "theiaviral"
  }
  if (! skip_rasusa) {
    # get genome length if it is not provided
    if (! defined(genome_length)) {
      call ncbi_datasets_task.ncbi_datasets_viral_taxon_summary as ncbi_taxon_summary {
        input:
          taxon = taxon
      }
    }
    # downsample reads to a specific coverage
    call rasusa_task.rasusa as rasusa {
      input:
        read1 = select_first([read_QC_trim.kraken2_extracted_read1]),
        read2 = select_first([read_QC_trim.kraken2_extracted_read2]),
        samplename = samplename,
        genome_length = select_first([genome_length, ncbi_taxon_summary.avg_genome_length])
    }
  }
  # clean read screening
  if (! skip_screen) {
    call read_screen_task.check_reads as clean_check_reads {
      input:
        read1 = select_first([rasusa.read1_subsampled, read_QC_trim.kraken2_extracted_read1]),
        read2 = select_first([rasusa.read2_subsampled, read_QC_trim.kraken2_extracted_read2]),
        workflow_series = "theiaviral"
    }
  }
  if (select_first([clean_check_reads.read_screen, ""]) == "PASS" || skip_screen) {
    # run de novo if no reference genome is provided so we can select a reference
    if (! defined(reference_fasta)) {
      # de novo assembly - prioritize metaviralspades
      if (! skip_metaviralspades) {
        call spades_task.spades {
          input:
            read1 = select_first([rasusa.read1_subsampled, read_QC_trim.kraken2_extracted_read1]),
            read2 = select_first([rasusa.read2_subsampled, read_QC_trim.kraken2_extracted_read2]),
            samplename = samplename,
            spades_type = "metaviral"
        }
      }
      # fallback to megahit if metaviralspades fails to identify a complete virus
      if (select_first([spades.spades_status, "FAIL"]) == "FAIL") {
        call megahit_task.megahit {
          input:
            read1 = select_first([rasusa.read1_subsampled, read_QC_trim.kraken2_extracted_read1]),
            read2 = select_first([rasusa.read2_subsampled, read_QC_trim.kraken2_extracted_read2]),
            samplename = samplename
        }
      }
      # quality control metrics for de novo assembly (ie. completeness, viral gene count, contamination)
      call checkv_task.checkv as checkv_denovo {
        input:
          assembly = select_first([spades.assembly_fasta, megahit.assembly_fasta]),
          samplename = samplename
      }
      # quality control metrics for de novo assembly (ie. contigs, n50, GC content, genome length)
      call quast_task.quast as quast_denovo {
        input:
          assembly = select_first([spades.assembly_fasta, megahit.assembly_fasta]),
          samplename = samplename
      }
      # ANI-based reference genome selection
      call skani_task.skani as skani {
        input:
          assembly_fasta = select_first([spades.assembly_fasta, megahit.assembly_fasta]),
          samplename = samplename
      }
      # download the best reference determined from skani
      call ncbi_datasets_task.ncbi_datasets_download_genome_accession as ncbi_datasets {
        input:
          ncbi_accession = skani.skani_top_accession,
          use_ncbi_virus = true
      }
    }
    # align reads to reference
    call bwa_task.bwa {
      input:
        samplename = samplename,
        read1 = select_first([rasusa.read1_subsampled, read_QC_trim.kraken2_extracted_read1]),
        read2 = select_first([rasusa.read2_subsampled, read_QC_trim.kraken2_extracted_read2]),
        reference_genome = select_first([reference_fasta, ncbi_datasets.ncbi_datasets_assembly_fasta])
    }
    # consensus calling via ivar
    call ivar_consensus.consensus {
      input:
        bamfile = bwa.sorted_bam,
        samplename = samplename,
        reference_genome = select_first([reference_fasta, ncbi_datasets.ncbi_datasets_assembly_fasta]),
        min_qual = min_map_quality,
        consensus_min_depth = select_first([min_depth, 10]),
        consensus_min_freq = min_allele_freq,
        all_positions = true
    }
    # variant calling via ivar
    call variant_call_task.variant_call as ivar_variants {
      input:
        mpileup = consensus.sample_mpileup,
        samplename = samplename,
        reference_genome = select_first([reference_fasta, ncbi_datasets.ncbi_datasets_assembly_fasta]),
        min_qual = min_map_quality,
        organism = "",
        variant_min_freq = min_allele_freq,
        variant_min_depth = select_first([min_depth, 10])
    }
    # quality control metrics for reads mapping to reference (ie. coverage, depth, base/map quality)
    call assembly_metrics_task.stats_n_coverage as read_mapping_stats {
      input:
        bamfile = bwa.sorted_bam,
        samplename = samplename
    }
    # quality control metrics for consensus (ie. number of bases, degenerate bases, genome length)
    call consensus_qc_task.consensus_qc as consensus_qc {
      input:
        assembly_fasta = consensus.consensus_seq,
        reference_genome = select_first([reference_fasta, ncbi_datasets.ncbi_datasets_assembly_fasta])
    }
    # quality control metrics for consensus (ie. completeness, viral gene count, contamination)
    call checkv_task.checkv as checkv_consensus {
      input:
        assembly = consensus.consensus_seq,
        samplename = samplename
    }
  }
  output {
    # versioning outputs
    String? theiaviral_illumina_pe_version = version_capture.phb_version
    String? theiaviral_illumina_pe_date = version_capture.date
    # ncbi datasets - taxon identification
    String ncbi_identify_taxon_id = ncbi_identify.taxon_id
    String ncbi_identify_taxon_name = ncbi_identify.taxon_name
    String ncbi_identify_read_extraction_rank = ncbi_identify.taxon_rank
    String ncbi_datasets_version = ncbi_identify.ncbi_datasets_version
    String ncbi_datasets_docker = ncbi_identify.ncbi_datasets_docker    
    # ncbi datasets - taxon summary
    File? ncbi_taxon_summary_tsv = ncbi_taxon_summary.taxon_summary_tsv
    Int? ncbi_taxon_summary_avg_genome_length = ncbi_taxon_summary.avg_genome_length
    # raw read quality control
    Int? fastq_scan_num_reads_raw1 = read_QC_trim.fastq_scan_raw1
    Int? fastq_scan_num_reads_raw2 = read_QC_trim.fastq_scan_raw2
    String? fastq_scan_raw_pairs = read_QC_trim.fastq_scan_raw_pairs
    String? fastq_scan_version = read_QC_trim.fastq_scan_version
    String? fastq_scan_docker = read_QC_trim.fastq_scan_docker
    File? fastq_scan_raw1_json = read_QC_trim.fastq_scan_raw1_json
    File? fastq_scan_raw2_json = read_QC_trim.fastq_scan_raw2_json
    # NCBI scrubbing
    File? read1_dehosted = read_QC_trim.read1_dehosted
    File? read2_dehosted = read_QC_trim.read2_dehosted
    Int? ncbi_scrub_human_spots_removed = read_QC_trim.ncbi_scrub_human_spots_removed
    String? ncbi_scrub_docker = read_QC_trim.ncbi_scrub_docker
    # trimming outputs - adapter trimming
    String? trimmomatic_version = read_QC_trim.trimmomatic_version
    String? trimmomatic_docker = read_QC_trim.trimmomatic_docker
    String? fastp_version = read_QC_trim.fastp_version
    File? fastp_html_report = read_QC_trim.fastp_html_report
    # bbduk outputs
    String? bbduk_docker = read_QC_trim.bbduk_docker
    File? bbduk_read1_clean = read_QC_trim.read1_clean
    File? bbduk_read2_clean = read_QC_trim.read2_clean
    # kraken2 outputs - taxonomic classification and read extraction
    String? kraken_version = read_QC_trim.kraken_version
    String? kraken_docker = read_QC_trim.kraken_docker
    String? kraken_database = read_QC_trim.kraken_database
    String? kraken_report = read_QC_trim.kraken_report
    File? kraken2_extracted_read1 = read_QC_trim.kraken2_extracted_read1
    File? kraken2_extracted_read2 = read_QC_trim.kraken2_extracted_read2
    # clean read quality control
    Int? fastq_scan_num_reads_clean1 = read_QC_trim.fastq_scan_clean1
    Int? fastq_scan_num_reads_clean2 = read_QC_trim.fastq_scan_clean2
    String? fastq_scan_clean_pairs = read_QC_trim.fastq_scan_clean_pairs
    File? fastq_scan_clean1_json = read_QC_trim.fastq_scan_clean1_json
    File? fastq_scan_clean2_json = read_QC_trim.fastq_scan_clean2_json
    # denovo genome assembly
    File? assembly_denovo_fasta = select_first([spades.assembly_fasta, megahit.assembly_fasta])
    String? metaviralspades_status = spades.spades_status
    String? metaviralspades_version = spades.spades_version
    String? metaviralspades_docker = spades.spades_docker
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
    String? quast_denovo_pipeline_date = quast_denovo.pipeline_date
    String? quast_denovo_version = quast_denovo.version
    String? quast_denovo_docker = quast_denovo.quast_docker
    # skani outputs - ANI-based reference genome selection
    File? skani_report = skani.skani_report
    String? skani_top_accession = skani.skani_top_accession
    Float? skani_top_score = skani.skani_top_score
    Float? skani_top_ani = skani.skani_top_ani
    Float? skani_top_ref_coverage = skani.skani_top_ref_coverage
    String? skani_database = skani.skani_database
    String? skani_version = skani.skani_version
    String? skani_docker = skani.skani_docker
    # ncbi_datasets outputs - download reference genome
    File? skani_top_ani_fasta = ncbi_datasets.ncbi_datasets_assembly_fasta
    # bwa outputs - reads aligned to best reference
    String? bwa_version = bwa.bwa_version
    String? bwa_samtools_version = bwa.sam_version
    File? bwa_read1_aligned = bwa.read1_aligned
    File? bwa_read2_aligned = bwa.read2_aligned
    File? bwa_sorted_bam = bwa.sorted_bam
    File? bwa_aligned_bai = bwa.sorted_bai
    File? bwa_read1_unaligned = bwa.read1_unaligned
    File? bwa_read2_unaligned = bwa.read2_unaligned
    File? sorted_bam_unaligned = bwa.sorted_bam_unaligned
    File? sorted_bam_unaligned_bai = bwa.sorted_bam_unaligned_bai
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
    Int consensus_n_variant_min_depth = select_first([min_depth, 20])
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
  }
}