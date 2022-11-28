version 1.0

import "../tasks/quality_control/task_trimmomatic.wdl" as trimmomatic
import "../tasks/quality_control/task_bbduk.wdl" as bbduk
import "../tasks/quality_control/task_fastq_scan.wdl" as fastq_scan
import "../tasks/taxon_id/task_midas.wdl" as midas

workflow read_QC_trim {
  meta {
    description: "Runs basic QC (fastq_scan), trimming (Trimmomatic), and adapter removal (bbduk) on illumina SE reads"
  }

  input {
    String samplename
    File read1_raw
    Int? trimmomatic_minlen = 25
    Int? trimmomatic_quality_trim_score = 30
    Int? trimmomatic_window_size = 4
    Int  bbduk_mem = 8
    Boolean call_midas = false
    File?    midas_db
  }
#  call read_clean.ncbi_scrub_se {
#    input:
#      samplename = samplename,
#      read1 = read1_raw
#  }
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
    call midas.midas as midas {
      input:
        samplename = samplename,
        read1 = read1_raw,
        midas_db = midas_db
    }
  }
#  call taxonID.kraken2 as kraken2_raw {
#    input:
#      samplename = samplename,
#      read1 = bbduk_se.read1_clean
#  }
#  call taxonID.kraken2 as kraken2_dehosted {
#    input:
#      samplename = samplename,
#      read1 = ncbi_scrub_se.read1_dehosted
#  }

  output {
    File read1_clean = bbduk_se.read1_clean

    Int fastq_scan_raw_number_reads = fastq_scan_raw.read1_seq
    Int fastq_scan_clean_number_reads = fastq_scan_clean.read1_seq

#    String  kraken_version            = kraken2_raw.version
#    Float   kraken_human              = kraken2_raw.percent_human
#    Float   kraken_sc2                = kraken2_raw.percent_sc2
#    String  kraken_report             = kraken2_raw.kraken_report
#    Float    kraken_human_dehosted    =    kraken2_dehosted.percent_human
#    Float    kraken_sc2_dehosted    =    kraken2_dehosted.percent_sc2
#    String    kraken_report_dehosted    =    kraken2_dehosted.kraken_report

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