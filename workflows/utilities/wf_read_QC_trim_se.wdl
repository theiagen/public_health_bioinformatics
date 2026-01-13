version 1.0

import "../../tasks/quality_control/basic_statistics/task_fastq_scan.wdl" as fastq_scan
import "../../tasks/quality_control/basic_statistics/task_fastqc.wdl" as fastqc_task
import "../../tasks/quality_control/read_filtering/task_bbduk.wdl" as bbduk_task
import "../../tasks/quality_control/read_filtering/task_fastp.wdl" as fastp_task
import "../../tasks/quality_control/read_filtering/task_trimmomatic.wdl" as trimmomatic
import "../../tasks/quality_control/read_filtering/task_ncbi_scrub.wdl" as ncbi_scrub
import "../../tasks/taxon_id/contamination/task_kraken2.wdl" as kraken
import "../../tasks/taxon_id/contamination/task_midas.wdl" as midas_task

workflow read_QC_trim_se {
  meta {
    description: "Runs basic QC (fastq-scan), trimming (SeqyClean), and taxonomic ID (Kraken2) on illumina SE reads"
  }
  input {
    String samplename
    File read1
    Int trim_min_length = 25
    Int trim_quality_min_score = 30
    Int trim_window_size = 4
    Int bbduk_memory = 8
    String? target_organism
    File? adapters
    File? phix
    String? workflow_series
    String? trimmomatic_args
    Boolean call_midas = false
    File? midas_db
    Boolean call_kraken = false
    File? kraken_db
    Int? kraken_disk_size
    Int? kraken_memory
    Int? kraken_cpu
    String read_processing = "trimmomatic" # options: trimmomatic, fastp
    String read_qc = "fastq_scan" # options: fastq_scan, fastqc
    String fastp_args = "-g -5 20 -3 20"
  }

  if (("~{workflow_series}" == "theiacov") || ("~{workflow_series}" == "theiameta")) {
    call ncbi_scrub.ncbi_scrub_se {
      input:
        samplename = samplename,
        read1 = read1
    }
  }

  if (read_processing == "trimmomatic") {
    call trimmomatic.trimmomatic_se {
      input:
        samplename = samplename,
        read1 = select_first([ncbi_scrub_se.read1_dehosted, read1]),
        trimmomatic_min_length = trim_min_length,
        trimmomatic_quality_trim_score = trim_quality_min_score,
        trimmomatic_window_size = trim_window_size,
        trimmomatic_args = trimmomatic_args
    }
  }
  if (read_processing == "fastp") {
    call fastp_task.fastp_se {
      input:
        samplename = samplename,
        read1 = select_first([ncbi_scrub_se.read1_dehosted, read1]),
        fastp_window_size = trim_window_size,
        fastp_quality_trim_score = trim_quality_min_score,
        fastp_min_length = trim_min_length,
        fastp_args = fastp_args
    }
  }
  call bbduk_task.bbduk_se {
    input:
      samplename = samplename,
      read1_trimmed = select_first([trimmomatic_se.read1_trimmed, fastp_se.read1_trimmed]),
      memory = bbduk_memory,
      adapters = adapters,
      phix = phix
  }
  if (read_qc == "fastq_scan") {
    call fastq_scan.fastq_scan as fastq_scan_raw {
      input:
        read1 = read1
    }
    call fastq_scan.fastq_scan as fastq_scan_clean {
      input:
        read1 = bbduk_se.read1_clean
    }
  }
  if (read_qc == "fastqc") {
    call fastqc_task.fastqc_se as fastqc_raw {
      input:
        read1 = read1
    }
    call fastqc_task.fastqc_se as fastqc_clean {
      input:
        read1 = bbduk_se.read1_clean
    }
  }
  if ("~{workflow_series}" == "theiacov") {
    call kraken.kraken2_theiacov as kraken2_theiacov_raw {
      input:
        samplename = samplename,
        read1 = read1,
        target_organism = target_organism,
        kraken2_db = kraken_db,
        disk_size = kraken_disk_size,
        memory = kraken_memory,
        cpu = kraken_cpu
    }
    call kraken.kraken2_theiacov as kraken2_theiacov_dehosted {
      input:
        samplename = samplename,
        read1 = select_first([ncbi_scrub_se.read1_dehosted]),
        target_organism = target_organism,
        kraken2_db = kraken_db,
        disk_size = kraken_disk_size,
        memory = kraken_memory,
        cpu = kraken_cpu
    }
  }
  if ("~{workflow_series}" == "theiaprok") {
    if (call_midas) {
      call midas_task.midas {
        input:
          samplename = samplename,
          read1 = read1,
          midas_db = midas_db
      }
    }
  }
  if ("~{workflow_series}" == "theiaprok") {
    if ((call_kraken) && defined(kraken_db)) {
      call kraken.kraken2_standalone {
        input:
          samplename = samplename,
          read1 = read1,
          kraken2_db = select_first([kraken_db]),
          disk_size = kraken_disk_size,
          memory = kraken_memory,
          cpu = kraken_cpu
      }
    }  if ((call_kraken) && ! defined(kraken_db)) {
      String kraken_db_warning = "Kraken database not defined"
    }
  }
  output {
    # NCBI scrubber
    File? read1_dehosted = ncbi_scrub_se.read1_dehosted
    Int? ncbi_scrub_human_spots_removed = ncbi_scrub_se.human_spots_removed
    String? ncbi_scrub_docker = ncbi_scrub_se.ncbi_scrub_docker

    # bbduk
    File read1_clean = bbduk_se.read1_clean
    String bbduk_docker = bbduk_se.bbduk_docker

    # fastq_scan raw (per read stats)
    Int? fastq_scan_raw1 = fastq_scan_raw.read1_seq
    Float? fastq_scan_r1_mean_readlength_raw = fastq_scan_raw.read1_mean_length
    Float? fastq_scan_r1_mean_q_raw = fastq_scan_raw.read1_mean_quality
    String? fastq_scan_version = fastq_scan_raw.version
    String? fastq_scan_docker = fastq_scan_raw.fastq_scan_docker
    File? fastq_scan_raw1_json = fastq_scan_raw.read1_fastq_scan_json

    # fastq_scan clean (per read stats)
    Int? fastq_scan_clean1 = fastq_scan_clean.read1_seq
    Float? fastq_scan_r1_mean_readlength_clean = fastq_scan_clean.read1_mean_length
    Float? fastq_scan_r1_mean_q_clean = fastq_scan_clean.read1_mean_quality
    File? fastq_scan_clean1_json = fastq_scan_clean.read1_fastq_scan_json

    # fastqc
    Int? fastqc_raw1 = fastqc_raw.read1_seq
    Int? fastqc_clean1 = fastqc_clean.read1_seq
    String? fastqc_version = fastqc_raw.version
    String? fastqc_docker = fastqc_raw.fastqc_docker
    File? fastqc_raw1_html = fastqc_raw.read1_fastqc_html
    File? fastqc_clean1_html = fastqc_clean.read1_fastqc_html
    
    # kraken2 - raw and dehosted
    String kraken_version = select_first([kraken2_theiacov_raw.version, kraken2_standalone.kraken2_version, ""])
    Float? kraken_human = kraken2_theiacov_raw.percent_human
    String? kraken_sc2 = kraken2_theiacov_raw.percent_sc2
    String? kraken_target_organism = kraken2_theiacov_raw.percent_target_organism
    String kraken_report = select_first([kraken2_theiacov_raw.kraken_report, kraken2_standalone.kraken2_report, ""])
    Float? kraken_human_dehosted = kraken2_theiacov_dehosted.percent_human
    String? kraken_sc2_dehosted = kraken2_theiacov_dehosted.percent_sc2
    String? kraken_target_organism_dehosted = kraken2_theiacov_dehosted.percent_target_organism
    String? kraken_target_organism_name = target_organism
    File? kraken_report_dehosted = kraken2_theiacov_dehosted.kraken_report
    String kraken_docker = select_first([kraken2_theiacov_raw.docker, kraken2_standalone.kraken2_docker, ""])
    String kraken_database = select_first([kraken2_theiacov_raw.database, kraken2_standalone.kraken2_database, kraken_db_warning, ""])
   
    # trimming versioning
    String? trimmomatic_version = trimmomatic_se.version
    String? trimmomatic_docker = trimmomatic_se.trimmomatic_docker
    String? fastp_version = fastp_se.version
    File? fastp_html_report = fastp_se.fastp_stats
    
    # midas
    String? midas_docker = midas.midas_docker
    File? midas_report = midas.midas_report
    String? midas_primary_genus = midas.midas_primary_genus
    String? midas_secondary_genus = midas.midas_secondary_genus
    Float? midas_secondary_genus_abundance = midas.midas_secondary_genus_abundance
    Float? midas_secondary_genus_coverage = midas.midas_secondary_genus_coverage
  }
}