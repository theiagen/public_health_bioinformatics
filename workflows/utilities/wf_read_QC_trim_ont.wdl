version 1.0

import "../../tasks/quality_control/task_fastq_scan.wdl" as fastq_scan
import "../../tasks/quality_control/task_nanoplot.wdl" as nanoplot_task
import "../../tasks/quality_control/task_nanoq.wdl" as nanoq_task
import "../../tasks/utilities/task_rasusa.wdl" as rasusa_task
import "../../tasks/utilities/task_kmc.wdl" as kmc_task
import "../../tasks/gene_typing/task_tiptoft.wdl" as tiptoft_task

workflow read_QC_trim_ont {
  meta {
    description: "Runs basic QC on Oxford Nanopore (ONT) reads with (1) fastq_scan, (2) nanoplot, (3) kmc, (4) rasusa downsampling, (5) tiptoft plasmid detection, and (6) nanoq filtering"
  }
  input {
    String samplename
    File reads
  }
  # perform fastq-scan on raw reads
  call fastq_scan.fastq_scan_se as fastq_scan_raw {
    input:
      read1 = reads,
      read1_name = samplename
  }
  # nanoplot
  call nanoplot_task.nanoplot {
    input:
      reads = reads,
      samplename = samplename
  }
  # kmc for genome size estimation
  call kmc_task.kmc {
    input:
      read1 = reads,
      samplename = samplename
  }
  # rasusa for random downsampling
  call rasusa_task.rasusa {
    input:
      read1 = reads,
      samplename = samplename,
      coverage = 150,
      genome_size = kmc.est_genome_size
  }
  # tiptoft
  call tiptoft_task.tiptoft {
    input:
      reads = reads,
      samplename = samplename
  }
  # nanoq/filtlong (default min length 500)
  call nanoq_task.nanoq {
    input:
      reads = rasusa.read1_subsampled,
      samplename = samplename
  }
  # perform fastq_scan again after cleaning for comparison
  call fastq_scan.fastq_scan_se as fastq_scan_clean {
    input:
      read1 = nanoq.filtered_reads,
      read1_name = samplename
  }
  output {
    # fastq scan raw outputs
    Int number_raw_reads = fastq_scan_raw.read1_seq

    # nanoplot outputs    
    File nanoplot_html = nanoplot.nanoplot_html
    String nanoplot_version = nanoplot.nanoplot_version
    
    # kmc outputs
    String est_genome_size = kmc.est_genome_size
    File kmc_kmer_stats = kmc.kmer_stats
    String kmc_version = kmc.kmc_version

    # rasusa outputs
    String rasusa_version = rasusa.rasusa_version  

    # tiptoft outputs
    File tiptoft_plasmid_replicon_fastq = tiptoft.tiptoft_plasmid_replicon_fastq
    File tiptoft_result_tsv = tiptoft.tiptoft_tsv
    String tiptoft_plasmid_replicon_genes = tiptoft.plasmid_replicon_genes
    String tiptoft_version = tiptoft.tiptoft_version

    # nanoq outputs
    File reads_clean = nanoq.filtered_reads
    String nanoq_version = nanoq.version

    # fastq scan clean outputs
    File fastq_scan_report = fastq_scan_clean.fastq_scan_report
    String fastq_scan_version = fastq_scan_clean.version
    Int number_clean_reads = fastq_scan_clean.read1_seq
  }
}