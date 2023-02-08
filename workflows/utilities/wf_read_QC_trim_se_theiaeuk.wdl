version 1.0

import "../../tasks/quality_control/task_trimmomatic.wdl" as trimmomatic
import "../../tasks/quality_control/task_bbduk.wdl" as bbduk
import "../../tasks/quality_control/task_fastq_scan.wdl" as fastq_scan
import "../../tasks/taxon_id/task_midas.wdl" as midas_task

workflow read_QC_trim_se {
  meta {
    description: "Runs basic QC (fastq_scan), trimming (Trimmomatic), and adapter removal (bbduk) on illumina SE reads"
  }

  input {
    String samplename
    File read1_raw
    Int trimmomatic_minlen = 25
    Int trimmomatic_quality_trim_score = 30
    Int trimmomatic_window_size = 4
    Int  bbduk_mem = 8
    Boolean call_midas = false
    File?    midas_db
  }
  call trimmomatic.trimmomatic_se {
    input:
      samplename = samplename,
      read1 = read1_raw,
      trimmomatic_minlen = trimmomatic_minlen,
      trimmomatic_quality_trim_score = trimmomatic_quality_trim_score,
      trimmomatic_window_size = trimmomatic_window_size
  }
  call bbduk.bbduk_se {
    input:
      samplename = samplename,
      read1_trimmed = trimmomatic_se.read1_trimmed,
      mem_size_gb = bbduk_mem
  }
  call fastq_scan.fastq_scan_se as fastq_scan_raw {
    input:
      read1 = read1_raw
  }
  call fastq_scan.fastq_scan_se as fastq_scan_clean {
    input:
      read1 = bbduk_se.read1_clean
  }
  if (call_midas) {
    call midas_task.midas {
      input:
        samplename = samplename,
        read1 = read1_raw,
        midas_db = midas_db
    }
  }
  output {
    File read1_clean = bbduk_se.read1_clean

    Int fastq_scan_raw_number_reads = fastq_scan_raw.read1_seq
    Int fastq_scan_clean_number_reads = fastq_scan_clean.read1_seq
    String fastq_scan_version = fastq_scan_raw.version
  
    String bbduk_docker = bbduk_se.bbduk_docker
    String trimmomatic_version = trimmomatic_se.version

    String? midas_docker = midas.midas_docker
    File? midas_report = midas.midas_report
    String? midas_primary_genus = midas.midas_primary_genus
    String? midas_secondary_genus = midas.midas_secondary_genus
    String? midas_secondary_genus_coverage = midas.midas_secondary_genus_coverage
  }
}