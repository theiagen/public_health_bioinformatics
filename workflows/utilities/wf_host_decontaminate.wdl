version 1.0

import "../../tasks/utilities/data_import/task_ncbi_datasets.wdl" as ncbi_datasets
import "../../tasks/alignment/task_minimap2.wdl" as minimap2_task
import "../../tasks/utilities/data_handling/task_parse_mapping.wdl" as parse_mapping_task
import "../../tasks/quality_control/basic_statistics/task_assembly_metrics.wdl" as assembly_metrics_task
import "../../tasks/taxon_id/task_identify_taxon_id.wdl" as identify_taxon_id_task

workflow host_decontaminate {
  meta {
    description: "Remove host reads from sequencing data via read mapping"
  }
  input {
    String samplename
    File read1
    File? read2
    String host
    Boolean is_accession = false
    Boolean refseq = true
    Boolean complete_only = false
    Int minimap2_memory = 32
  }
  String hostsample = samplename + "_host"
  # gather an accession from a taxon input
  if (! is_accession) {
    call identify_taxon_id_task.identify_taxon_id as ncbi_identify {
      input:
        taxon = host,
        refseq = refseq,
        complete = complete_only,
        summary_limit = 1,
        use_ncbi_virus = false
    }
  }
  # download accession
  call ncbi_datasets.ncbi_datasets_download_genome_accession as download_accession {
    input:
      ncbi_accession = select_first([ncbi_identify.ncbi_datasets_accession, host]),
      use_ncbi_virus = false
  }
  # run paired-end mode
  if (defined(read2)) {
    call minimap2_task.minimap2 as minimap2_pe {
      input:
        samplename = hostsample,
        query1 = read1,
        query2 = read2,
        reference = select_first([download_accession.ncbi_datasets_assembly_fasta]),
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
        samplename = hostsample,
        query1 = read1,
        reference = select_first([download_accession.ncbi_datasets_assembly_fasta]),
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
      samplename = hostsample
  }
  # extract unaligned reads from BAM (those that don't map to host)
  call parse_mapping_task.bam_to_unaligned_fastq {
    input:
      bam = parse_mapping.bam,
      samplename = hostsample,
      paired = defined(read2)
  }
  # calculate read mapping statistics
  call assembly_metrics_task.stats_n_coverage as read_mapping_stats {
    input:
      bamfile = parse_mapping.bam,
      samplename = hostsample
  }
  output {
    # Datasets download outputs
    File? host_genome_fasta = download_accession.ncbi_datasets_assembly_fasta
    File? host_genome_data_report_json = download_accession.ncbi_datasets_assembly_data_report_json
    String host_genome_accession = select_first([ncbi_identify.ncbi_datasets_accession, host])
    String ncbi_datasets_version = download_accession.ncbi_datasets_version
    # Read mapping outputs
    File? host_mapped_sorted_bam = parse_mapping.bam
    File? host_mapped_sorted_bai = parse_mapping.bai
    File? dehost_read1 = bam_to_unaligned_fastq.read1_unaligned
    File? dehost_read2 = bam_to_unaligned_fastq.read2_unaligned
    File? host_bam = parse_mapping.bam
    String? samtools_version = bam_to_unaligned_fastq.sam_version
    # Read mapping stats
    File? host_mapping_stats = read_mapping_stats.stats
    File? host_mapping_cov_hist = read_mapping_stats.cov_hist
    File? host_flagstat = read_mapping_stats.flagstat
    Float? host_mapping_coverage = read_mapping_stats.coverage
    Float? host_mapping_mean_depth = read_mapping_stats.depth
    Float? host_percent_mapped_reads = read_mapping_stats.percentage_mapped_reads
    File? host_mapping_metrics = read_mapping_stats.metrics_txt
  }
}