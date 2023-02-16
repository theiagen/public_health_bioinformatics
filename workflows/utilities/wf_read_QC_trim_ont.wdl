version 1.0

import "../../tasks/quality_control/task_fastq_scan.wdl" as fastq_scan
import "../../tasks/utilities/task_rasusa.wdl" as rasusa_task

workflow read_QC_trim_ont {
  meta {
    description: "Runs basic QC on Oxford Nanopore (ONT) reads"
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
  # call nanoplot_task.nanplot {

  # }

  # rasusa for random downsampling
  call rasusa_task.rasusa { # placeholder task
    input:
      read1 = reads,
      samplename = samplename,
      coverage = 150,
      genome_size = 5000000 # placeholder values
  }

  # tiptoft
  # call tiptoft_task.tiptoft {

  # }
  # porechop? (optional)
  # call porechop_task.porechop {

  # }
  # nanoq/filtlong (default min length 500)
  # call nanoq_task.nanoq {

  # }
 
  # perform fastq_scan again after cleaning for comparison
  call fastq_scan.fastq_scan_se as fastq_scan_clean {
    input:
      read1 = reads,
      read1_name = samplename
  }
  output {
    # placeholder outputs
    File reads_clean = reads # we're not cleaning anything, so return the input reads
    Int number_raw_reads = fastq_scan_raw.read1_seq
    File fastq_scan_report = fastq_scan_clean.fastq_scan_report
    String fastq_scan_version = fastq_scan_clean.version
    Int number_clean_reads = fastq_scan_clean.read1_seq
    String rasusa_version = rasusa.rasusa_version
    
  }
}