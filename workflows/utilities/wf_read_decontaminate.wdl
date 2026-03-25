version 1.0

import "../../tasks/utilities/data_import/task_ncbi_datasets.wdl" as ncbi_datasets
import "../../tasks/alignment/task_minimap2.wdl" as minimap2_task
import "../../tasks/utilities/data_handling/task_parse_mapping.wdl" as parse_mapping_task
import "../../tasks/quality_control/basic_statistics/task_mapping_stats.wdl" as mapping_stats_task
import "../../tasks/utilities/task_datasets_genome_length.wdl" as identify_genome_task
import "../../tasks/taxon_id/contamination/task_contaminant_check.wdl" as contaminant_check_task

workflow read_decontaminate {
  meta {
    description: "Remove contaminant reads from sequencing data via read mapping"
  }
  input {
    String samplename
    File read1
    File? read2
    String contaminant
    Boolean is_accession = false
    Boolean is_genome = false
    Boolean refseq = true
    Boolean complete_only = false

    String? expected_sequences
    Float min_expected_coverage = 0
    Int min_expected_depth = 0
    Int minimap2_memory = 32
  }
  String contaminant_samplename = samplename + "_contaminant"
  # gather an accession from a taxon input
  if (! is_accession && ! is_genome) {
    call identify_genome_task.datasets_genome_length as ncbi_identify {
      input:
        taxon = contaminant,
        refseq = refseq,
        complete = complete_only,
        summary_limit = 1,
        use_ncbi_virus = false
    }
  }
  # download accession if genome isn't directly provided
  if (! is_genome) {
    call ncbi_datasets.ncbi_datasets_download_genome_accession as download_accession {
      input:
        ncbi_accession = select_first([ncbi_identify.ncbi_datasets_accession, contaminant]),
        use_ncbi_virus = false
    }
  }
  # run paired-end mode
  if (defined(read2)) {
    call minimap2_task.minimap2 as minimap2_pe {
      input:
        samplename = contaminant_samplename,
        query1 = read1,
        query2 = read2,
        reference = select_first([download_accession.ncbi_datasets_assembly_fasta, contaminant]),
        mode = "sr",
        output_sam = true,
        long_read_flags = false,
        memory = minimap2_memory
    }
  }
  # run ONT mode (Illumina single-end reads are not supported)
  if (! defined(read2)) {
    call minimap2_task.minimap2 as minimap2_ont {
      input:
        samplename = contaminant_samplename,
        query1 = read1,
        reference = select_first([download_accession.ncbi_datasets_assembly_fasta, contaminant]),
        mode = "map-ont",
        output_sam = true,
        long_read_flags = true,
        memory = minimap2_memory
    }
  }
  # extract BAM
  call parse_mapping_task.sam_to_sorted_bam as parse_mapping {
    input:
      sam = select_first([minimap2_pe.minimap2_out, minimap2_ont.minimap2_out]),
      samplename = contaminant_samplename
  }
  # extract unaligned reads from BAM (those that don't map to host)
  call parse_mapping_task.bam_to_unaligned_fastq {
    input:
      bam = parse_mapping.bam,
      samplename = contaminant_samplename,
      paired = defined(read2)
  }
  # calculate read mapping statistics
  call mapping_stats_task.mapping_stats as read_mapping_stats {
    input:
      bamfile = parse_mapping.bam,
      samplename = contaminant_samplename,
      read1 = read1,
      read2 = read2
  }
  # run contaminant check
  if (defined(expected_sequences)) {
    if (read_mapping_stats.mapping_stats_status == "PASS") {
      call contaminant_check_task.contaminant_check {
        input:
          expected_sequences = select_first([expected_sequences]),
          coverage_by_sequence_json = select_first([read_mapping_stats.coverage_by_sequence_json]),
          depth_by_sequence_json = select_first([read_mapping_stats.depth_by_sequence_json]),
          min_percent_coverage = min_expected_coverage,
          min_depth = min_expected_depth,
          cov_stats = read_mapping_stats.cov_stats
      }
    }
    if (read_mapping_stats.mapping_stats_status == "FAIL") {
      String contaminant_check_fail = "FAIL: no reads mapped to inputted sequences"
      Map[String, Float] empty_cov = {"": 0}
      Map[String, Float] empty_depth = {"": 0}
    }
  }
  output {
    # Datasets download outputs
    File? contaminant_genome_fasta = download_accession.ncbi_datasets_assembly_fasta
    File? contaminant_genome_data_report_json = download_accession.ncbi_datasets_assembly_data_report_json
    String? contaminant_genome_accession = select_first([ncbi_identify.ncbi_datasets_accession, contaminant])
    String? ncbi_datasets_version = download_accession.ncbi_datasets_version
    # Read mapping outputs
    File? contaminant_mapped_sorted_bam = parse_mapping.bam
    File? contaminant_mapped_sorted_bai = parse_mapping.bai
    File? decontaminate_read1 = bam_to_unaligned_fastq.read1_unaligned
    File? decontaminate_read2 = bam_to_unaligned_fastq.read2_unaligned
    File? contaminant_bam = parse_mapping.bam
    String? samtools_version = bam_to_unaligned_fastq.sam_version
    # Read mapping stats
    File? contaminant_mapping_stats = read_mapping_stats.stats
    File? contaminant_mapping_cov_hist = read_mapping_stats.cov_hist
    File? contaminant_flagstat = read_mapping_stats.flagstat
    Float? contaminant_mapping_coverage = read_mapping_stats.coverage
    Float? contaminant_mapping_mean_depth = read_mapping_stats.depth
    Float? contaminant_percent_mapped_reads = read_mapping_stats.percentage_mapped_reads
    Map[String, Float]? contaminant_coverage_by_sequence = select_first([empty_cov, read_mapping_stats.coverage_by_sequence])
    Map[String, Float]? contaminant_depth_by_sequence = select_first([empty_depth, read_mapping_stats.depth_by_sequence])
    # Contaminant check outputs
    String? contaminant_check_status = select_first([contaminant_check_fail, contaminant_check.contaminant_check_status])
  }
}