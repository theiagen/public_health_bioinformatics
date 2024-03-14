version 1.0

import "../../tasks/quality_control/basic_statistics/task_fastq_scan.wdl" as fastq_scan
import "../../tasks/quality_control/basic_statistics/task_fastqc.wdl" as fastqc_task
import "../../tasks/quality_control/basic_statistics/task_readlength.wdl" as readlength_task
import "../../tasks/quality_control/read_filtering/task_bbduk.wdl" as bbduk_task
import "../../tasks/quality_control/read_filtering/task_fastp.wdl" as fastp_task
import "../../tasks/quality_control/read_filtering/task_ncbi_scrub.wdl" as ncbi_scrub
import "../../tasks/quality_control/read_filtering/task_trimmomatic.wdl" as trimmomatic
import "../../tasks/taxon_id/contamination/task_kraken2.wdl" as kraken
import "../../tasks/taxon_id/contamination/task_midas.wdl" as midas_task

workflow read_QC_trim_pe {
  meta {
    description: "Runs basic QC (fastq-scan), trimming (trimmomatic), and taxonomic ID (Kraken2) on illumina PE reads"
  }
  input {
    String samplename
    File read1
    File read2
    Int trim_min_length = 75
    Int trim_quality_min_score = 30
    Int trim_window_size = 4
    Int bbduk_memory = 8
    Boolean call_midas = true
    File? midas_db
    Boolean call_kraken = false
    Int? kraken_disk_size
    Int? kraken_memory
    File? kraken_db
    String? target_organism
    File? adapters
    File? phix
    String? workflow_series
    String read_processing = "trimmomatic" # options: trimmomatic, fastp
    String read_qc = "fastq_scan" # options: fastq_scan, fastqc
    String? trimmomatic_args
    String fastp_args = "--detect_adapter_for_pe -g -5 20 -3 20"
  }
  if (("~{workflow_series}" == "theiacov") || ("~{workflow_series}" == "theiameta")) {
    call ncbi_scrub.ncbi_scrub_pe {
      input:
        samplename = samplename,
        read1 = read1,
        read2 = read2
    }
  }
  if ("~{workflow_series}" == "theiacov") {
    call kraken.kraken2_theiacov as kraken2_theiacov_raw {
      input:
        samplename = samplename,
        read1 = read1,
        read2 = read2,
        target_organism = target_organism
    }
    call kraken.kraken2_theiacov as kraken2_theiacov_dehosted {
      input:
        samplename = samplename,
        read1 = select_first([ncbi_scrub_pe.read1_dehosted]),
        read2 = ncbi_scrub_pe.read2_dehosted,
        target_organism = target_organism
    }
  }
  if (read_processing == "trimmomatic"){
    call trimmomatic.trimmomatic_pe {
      input:
        samplename = samplename,
        read1 = select_first([ncbi_scrub_pe.read1_dehosted, read1]),
        read2 = select_first([ncbi_scrub_pe.read2_dehosted, read2]),
        trimmomatic_window_size = trim_window_size,
        trimmomatic_quality_trim_score = trim_quality_min_score,
        trimmomatic_min_length = trim_min_length,
        trimmomatic_args = trimmomatic_args
    }
  }
  if (read_processing == "fastp"){
    call fastp_task.fastp_pe as fastp {
      input:
        samplename = samplename,
        read1 = select_first([ncbi_scrub_pe.read1_dehosted, read1]),
        read2 = select_first([ncbi_scrub_pe.read2_dehosted, read2]),
        fastp_window_size = trim_window_size,
        fastp_quality_trim_score = trim_quality_min_score,
        fastp_min_length = trim_min_length,
        fastp_args = fastp_args
    }
  }
  call bbduk_task.bbduk {
    input:
      samplename = samplename,
      read1_trimmed = select_first([trimmomatic_pe.read1_trimmed, fastp.read1_trimmed]),
      read2_trimmed = select_first([trimmomatic_pe.read2_trimmed, fastp.read2_trimmed]),
      memory = bbduk_memory,
      adapters = adapters,
      phix = phix
  }
  if (read_qc == "fastqc") {
    call fastqc_task.fastqc as fastqc_raw {
      input:
        read1 = read1,
        read2 = read2
    }
    call fastqc_task.fastqc as fastqc_clean {
      input:
        read1 = bbduk.read1_clean,
        read2 = bbduk.read2_clean
    }
  }
  if (read_qc == "fastq_scan") {
    call fastq_scan.fastq_scan_pe as fastq_scan_raw {
      input:
        read1 = read1,
        read2 = read2,
    }
    call fastq_scan.fastq_scan_pe as fastq_scan_clean {
      input:
        read1 = bbduk.read1_clean,
        read2 = bbduk.read2_clean
    }
  }
  if ("~{workflow_series}" == "theiaprok") {
    if (call_midas) {
      call midas_task.midas {
        input:
          samplename = samplename,
          read1 = read1,
          read2 = read2,
          midas_db = midas_db
      }
    }
  }
  if ("~{workflow_series}" == "theiaprok") {
    if (call_kraken) {
      call kraken.kraken2_standalone {
        input:
          samplename = samplename,
          read1 = read1,
          read2 = read2,
          kraken2_db = select_first([kraken_db]),
          disk_size = kraken_disk_size,
          memory = kraken_memory
      }
    }
  }
  if ("~{workflow_series}" == "theiameta") {
    call readlength_task.readlength {
      input:
        read1 = bbduk.read1_clean,
        read2 = bbduk.read2_clean
    }
  }
  output {
    # NCBI scrubber
    File? read1_dehosted = ncbi_scrub_pe.read1_dehosted
    File? read2_dehosted = ncbi_scrub_pe.read2_dehosted
    Int? read1_human_spots_removed = ncbi_scrub_pe.read1_human_spots_removed
    Int? read2_human_spots_removed = ncbi_scrub_pe.read2_human_spots_removed
    String? ncbi_scrub_docker = ncbi_scrub_pe.ncbi_scrub_docker

    # bbduk
    File read1_clean = bbduk.read1_clean
    File read2_clean = bbduk.read2_clean
    String bbduk_docker = bbduk.bbduk_docker

    # fastq_scan
    Int? fastq_scan_raw1 = fastq_scan_raw.read1_seq
    Int? fastq_scan_raw2 = fastq_scan_raw.read2_seq
    String? fastq_scan_raw_pairs = fastq_scan_raw.read_pairs
    Int? fastq_scan_clean1 = fastq_scan_clean.read1_seq
    Int? fastq_scan_clean2 = fastq_scan_clean.read2_seq
    String? fastq_scan_clean_pairs = fastq_scan_clean.read_pairs
    String? fastq_scan_version = fastq_scan_raw.version
    String? fastq_scan_docker = fastq_scan_raw.fastq_scan_docker
    
    # fastqc
    Int? fastqc_raw1 = fastqc_raw.read1_seq
    Int? fastqc_raw2 = fastqc_raw.read2_seq
    String? fastqc_raw_pairs = fastqc_raw.read_pairs
    Int? fastqc_clean1 = fastqc_clean.read1_seq
    Int? fastqc_clean2 = fastqc_clean.read2_seq
    String? fastqc_clean_pairs = fastqc_clean.read_pairs
    String? fastqc_version = fastqc_raw.version
    String? fastqc_docker = fastqc_raw.fastqc_docker
    File? fastqc_raw1_html = fastqc_raw.read1_fastqc_html
    File? fastqc_raw2_html = fastqc_raw.read2_fastqc_html
    File? fastqc_clean1_html = fastqc_clean.read1_fastqc_html
    File? fastqc_clean2_html = fastqc_clean.read2_fastqc_html
    
    # kraken2 - theiacov and theiaprok
    String kraken_version = select_first([kraken2_theiacov_raw.version, kraken2_standalone.kraken2_version, ""])
    Float? kraken_human =  kraken2_theiacov_raw.percent_human
    Float? kraken_sc2 = kraken2_theiacov_raw.percent_sc2
    String? kraken_target_organism = kraken2_theiacov_raw.percent_target_organism
    String kraken_report = select_first([kraken2_theiacov_raw.kraken_report, kraken2_standalone.kraken2_report, ""])
    Float? kraken_human_dehosted = kraken2_theiacov_dehosted.percent_human
    Float? kraken_sc2_dehosted = kraken2_theiacov_dehosted.percent_sc2
    String? kraken_target_organism_dehosted = kraken2_theiacov_dehosted.percent_target_organism
    String? kraken_target_organism_name = target_organism
    File? kraken_report_dehosted = kraken2_theiacov_dehosted.kraken_report
    String kraken_docker = select_first([kraken2_theiacov_raw.docker, kraken2_standalone.kraken2_docker, ""])
    
    # trimming versioning
    String? trimmomatic_version = trimmomatic_pe.version
    String? trimmomatic_docker = trimmomatic_pe.trimmomatic_docker
    String? fastp_version = fastp.version
    File? fastp_html_report = fastp.fastp_stats

    # midas
    String? midas_docker = midas.midas_docker
    File? midas_report = midas.midas_report
    String? midas_primary_genus = midas.midas_primary_genus
    String? midas_secondary_genus = midas.midas_secondary_genus
    Float? midas_secondary_genus_abundance = midas.midas_secondary_genus_abundance
    Float? midas_secondary_genus_coverage = midas.midas_secondary_genus_coverage

    # readlength
    Float? average_read_length = readlength.average_read_length
  }
}