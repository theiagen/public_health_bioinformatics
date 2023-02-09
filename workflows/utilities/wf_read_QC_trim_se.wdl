version 1.0

import "../../tasks/quality_control/task_fastq_scan.wdl" as fastq_scan
import "../../tasks/quality_control/task_trimmomatic.wdl" as trimmomatic
import "../../tasks/quality_control/task_bbduk.wdl" as bbduk_task
import "../../tasks/quality_control/task_fastp.wdl" as fastp_task
import "../../tasks/taxon_id/task_kraken2.wdl" as kraken
import "../../tasks/taxon_id/task_midas.wdl" as midas_task

workflow read_QC_trim_se {
  meta {
    description: "Runs basic QC (fastq-scan), trimming (SeqyClean), and taxonomic ID (Kraken2) on illumina SE reads"
  }
  input {
    String samplename
    File read1_raw
    Int trim_minlen = 25
    Int trim_quality_trim_score = 30
    Int trim_window_size = 4
    Int bbduk_mem = 8
    String? target_org
    File? adapters
    File? phix
    String? workflow_series
    String? trimmomatic_args
    Boolean call_midas = false
    File? midas_db
    String read_processing = "trimmomatic"
    String fastp_args = "-g -5 20 -3 20"
  }
  if (read_processing == "trimmomatic"){
    call trimmomatic.trimmomatic_se {
      input:
        samplename = samplename,
        read1 = read1_raw,
        trimmomatic_minlen = trim_minlen,
        trimmomatic_quality_trim_score = trim_quality_trim_score,
        trimmomatic_window_size = trim_window_size,
        trimmomatic_args = trimmomatic_args
    }
  }
  if (read_processing == "fastp"){
    call fastp_task.fastp_se_theiaprok as fastp_se {
      input:
        samplename = samplename,
        read1 = read1_raw,
        fastp_window_size = trim_window_size,
        fastp_quality_trim_score = trim_quality_trim_score,
        fastp_minlen = trim_minlen,
        fastp_args = fastp_args
    }
  }
  call bbduk_task.bbduk_se {
    input:
      samplename = samplename,
      read1_trimmed = select_first([trimmomatic_se.read1_trimmed, fastp_se.read1_trimmed]),
      memory = bbduk_mem,
      adapters = adapters,
      phix = phix
  }
  call fastq_scan.fastq_scan_se as fastq_scan_raw {
    input:
      read1 = read1_raw
  }
  call fastq_scan.fastq_scan_se as fastq_scan_clean {
    input:
      read1 = bbduk_se.read1_clean
  }
  if ("~{workflow_series}" == "theiacov"){
    call kraken.kraken2_theiacov as kraken2_raw {
      input:
        samplename = samplename,
        read1 = bbduk_se.read1_clean,
        target_org = target_org
    }
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
    # bbduk
    File read1_clean = bbduk_se.read1_clean
    String bbduk_docker = bbduk_se.bbduk_docker

    # fastq_scan
    Int fastq_scan_raw_number_reads = fastq_scan_raw.read1_seq
    Int fastq_scan_clean_number_reads = fastq_scan_clean.read1_seq
    String fastq_scan_version = fastq_scan_raw.version
    
    # kraken2
    String? kraken_version = kraken2_raw.version
    Float? kraken_human = kraken2_raw.percent_human
    Float? kraken_sc2 = kraken2_raw.percent_sc2
    String? kraken_target_org = kraken2_raw.percent_target_org
    File? kraken_report = kraken2_raw.kraken_report
    String? kraken_target_org_name = target_org
   
    # trimming versioning
    String? trimmomatic_version = trimmomatic_se.version
    String? fastp_version = fastp_se.version
    
    # midas
    String? midas_docker = midas.midas_docker
    File? midas_report = midas.midas_report
    String? midas_primary_genus = midas.midas_primary_genus
    String? midas_secondary_genus = midas.midas_secondary_genus
    Float? midas_secondary_genus_abundance = midas.midas_secondary_genus_abundance
  }
}