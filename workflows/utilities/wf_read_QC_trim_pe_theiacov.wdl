version 1.0

import "../../tasks/quality_control/task_fastq_scan.wdl" as fastq_scan
import "../../tasks/quality_control/task_trimmomatic.wdl" as trimmomatic
import "../../tasks/quality_control/task_ncbi_scrub.wdl" as ncbi_scrub
import "../../tasks/quality_control/task_bbduk.wdl" as bbduk_task
import "../../tasks/taxon_id/task_kraken2.wdl" as kraken

workflow read_QC_trim_pe {
  meta {
    description: "Runs basic QC (fastq-scan), trimming (trimmomatic), and taxonomic ID (Kraken2) on illumina PE reads"
  }
  input {
    String samplename
    File read1_raw
    File read2_raw
    Int trimmomatic_minlen = 75
    Int trimmomatic_quality_trim_score = 30
    Int trimmomatic_window_size = 4
    Int bbduk_mem = 8
    String? target_org
    File? adapters
    File? phix
    String? trim_args
  }
  call ncbi_scrub.ncbi_scrub_pe {
    input:
      samplename = samplename,
      read1 = read1_raw,
      read2 = read2_raw
  }
  call trimmomatic.trimmomatic_pe {
    input:
      samplename = samplename,
      read1 = ncbi_scrub_pe.read1_dehosted,
      read2 = ncbi_scrub_pe.read2_dehosted,
      trimmomatic_minlen = trimmomatic_minlen,
      trimmomatic_quality_trim_score = trimmomatic_quality_trim_score,
      trimmomatic_window_size = trimmomatic_window_size,
      trimmomatic_args = trim_args
  }
  call bbduk_task.bbduk {
    input:
      samplename = samplename,
      read1_trimmed = trimmomatic_pe.read1_trimmed,
      read2_trimmed = trimmomatic_pe.read2_trimmed,
      memory = bbduk_mem,
      adapters = adapters,
      phix = phix
  }
  call fastq_scan.fastq_scan_pe as fastq_scan_raw {
    input:
      read1 = read1_raw,
      read2 = read2_raw,
  }
  call fastq_scan.fastq_scan_pe as fastq_scan_clean {
    input:
      read1 = bbduk.read1_clean,
      read2 = bbduk.read2_clean
  }
  call kraken.kraken2_theiacov as kraken2_raw {
    input:
      samplename = samplename,
      read1 = read1_raw,
      read2 = read2_raw,
      target_org = target_org
  }
  call kraken.kraken2_theiacov as kraken2_dehosted {
    input:
      samplename = samplename,
      read1 = ncbi_scrub_pe.read1_dehosted,
      read2 = ncbi_scrub_pe.read2_dehosted,
      target_org = target_org
  }
  output {
    File read1_dehosted = ncbi_scrub_pe.read1_dehosted
    File read2_dehosted = ncbi_scrub_pe.read2_dehosted
    Int read1_human_spots_removed = ncbi_scrub_pe.read1_human_spots_removed
    Int read2_human_spots_removed = ncbi_scrub_pe.read2_human_spots_removed
    File read1_clean = bbduk.read1_clean
    File read2_clean = bbduk.read2_clean
    Int fastq_scan_raw1 = fastq_scan_raw.read1_seq
    Int fastq_scan_raw2 = fastq_scan_raw.read2_seq
    String fastq_scan_raw_pairs = fastq_scan_raw.read_pairs
    Int fastq_scan_clean1 = fastq_scan_clean.read1_seq
    Int fastq_scan_clean2 = fastq_scan_clean.read2_seq
    String fastq_scan_clean_pairs = fastq_scan_clean.read_pairs
    String kraken_version = kraken2_raw.version
    Float kraken_human = kraken2_raw.percent_human
    Float kraken_sc2 = kraken2_raw.percent_sc2
    String? kraken_target_org = kraken2_raw.percent_target_org
    File kraken_report = kraken2_raw.kraken_report
    Float kraken_human_dehosted = kraken2_dehosted.percent_human
    Float kraken_sc2_dehosted = kraken2_dehosted.percent_sc2
    String? kraken_target_org_dehosted = kraken2_dehosted.percent_target_org
    String? kraken_target_org_name = target_org
    File kraken_report_dehosted = kraken2_dehosted.kraken_report
    String fastq_scan_version = fastq_scan_raw.version
    String bbduk_docker = bbduk.bbduk_docker
    String trimmomatic_version = trimmomatic_pe.version
  }
}