version 1.0

import "../../tasks/quality_control/task_fastq_scan.wdl" as fastq_scan
import "../../tasks/quality_control/task_fastqc.wdl" as fastqc_task
import "../../tasks/quality_control/task_trimmomatic.wdl" as trimmomatic
import "../../tasks/quality_control/task_ncbi_scrub.wdl" as ncbi_scrub
import "../../tasks/quality_control/task_bbduk.wdl" as bbduk_task
import "../../tasks/quality_control/task_readlength.wdl" as readlength_task
import "../../tasks/quality_control/task_fastp.wdl" as fastp_task
import "../../tasks/taxon_id/task_kraken2.wdl" as kraken
import "../../tasks/taxon_id/task_midas.wdl" as midas_task

workflow read_QC_trim_pe {
  meta {
    description: "Runs basic QC (fastq-scan), trimming (trimmomatic), and taxonomic ID (Kraken2) on illumina PE reads"
  }
  input {
    String samplename
    File read1_raw
    File read2_raw
    Int trim_minlen = 75
    Int trim_quality_trim_score = 30
    Int trim_window_size = 4
    Int bbduk_mem = 8
    Boolean call_midas = false
    File? midas_db
    Boolean call_kraken = false
    File? kraken_db
    String? target_org
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
        read1 = read1_raw,
        read2 = read2_raw
    }
  }
  if ("~{workflow_series}" == "theiacov") {
    call kraken.kraken2_theiacov as kraken2_theiacov_raw {
      input:
        samplename = samplename,
        read1 = read1_raw,
        read2 = read2_raw,
        target_org = target_org
    }
    call kraken.kraken2_theiacov as kraken2_theiacov_dehosted {
      input:
        samplename = samplename,
        read1 = select_first([ncbi_scrub_pe.read1_dehosted]),
        read2 = ncbi_scrub_pe.read2_dehosted,
        target_org = target_org
    }
  }
  if (read_processing == "trimmomatic"){
    call trimmomatic.trimmomatic_pe {
      input:
        samplename = samplename,
        read1 = select_first([ncbi_scrub_pe.read1_dehosted, read1_raw]),
        read2 = select_first([ncbi_scrub_pe.read2_dehosted, read2_raw]),
        trimmomatic_window_size = trim_window_size,
        trimmomatic_quality_trim_score = trim_quality_trim_score,
        trimmomatic_minlen = trim_minlen,
        trimmomatic_args = trimmomatic_args
    }
  }
  if (read_processing == "fastp"){
    call fastp_task.fastp_pe as fastp {
      input:
        samplename = samplename,
        read1 = select_first([ncbi_scrub_pe.read1_dehosted, read1_raw]),
        read2 = select_first([ncbi_scrub_pe.read2_dehosted, read2_raw]),
        fastp_window_size = trim_window_size,
        fastp_quality_trim_score = trim_quality_trim_score,
        fastp_minlen = trim_minlen,
        fastp_args = fastp_args
    }
  }
  call bbduk_task.bbduk {
    input:
      samplename = samplename,
      read1_trimmed = select_first([trimmomatic_pe.read1_trimmed, fastp.read1_trimmed]),
      read2_trimmed = select_first([trimmomatic_pe.read2_trimmed, fastp.read2_trimmed]),
      memory = bbduk_mem,
      adapters = adapters,
      phix = phix
  }
  if (read_qc == "fastqc") {
    call fastqc_task.fastqc as fastqc_raw {
      input:
        read1 = read1_raw,
        read2 = read2_raw
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
        read1 = read1_raw,
        read2 = read2_raw,
    }
    call fastq_scan.fastq_scan_pe as fastq_scan_clean {
      input:
        read1 = bbduk.read1_clean,
        read2 = bbduk.read2_clean
    }
  }
  if (call_midas) {
    call midas_task.midas {
      input:
        samplename = samplename,
        read1 = read1_raw,
        read2 = read2_raw,
        midas_db = midas_db
    }
  }
  if ("~{workflow_series}" == "theiaprok") {
    if (call_kraken) {
      call kraken.kraken2_standalone {
        input:
          samplename = samplename,
          read1 = read1_raw,
          read2 = read2_raw,
          kraken2_db = select_first([kraken_db])
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
    String? kraken_target_org = kraken2_theiacov_raw.percent_target_org
    String kraken_report = select_first([kraken2_theiacov_raw.kraken_report, kraken2_standalone.kraken2_report, ""])
    Float? kraken_human_dehosted = kraken2_theiacov_dehosted.percent_human
    Float? kraken_sc2_dehosted = kraken2_theiacov_dehosted.percent_sc2
    String? kraken_target_org_dehosted = kraken2_theiacov_dehosted.percent_target_org
    String? kraken_target_org_name = target_org
    File? kraken_report_dehosted = kraken2_theiacov_dehosted.kraken_report
    String kraken_docker = select_first([kraken2_theiacov_raw.docker, kraken2_standalone.kraken2_docker, ""])
    
    # trimming versioning
    String? trimmomatic_version = trimmomatic_pe.version
    String? trimmomatic_docker = trimmomatic_pe.trimmomatic_docker
    String? fastp_version = fastp.version

    # midas
    String? midas_docker = midas.midas_docker
    File? midas_report = midas.midas_report
    String? midas_primary_genus = midas.midas_primary_genus
    String? midas_secondary_genus = midas.midas_secondary_genus
    Float? midas_secondary_genus_abundance = midas.midas_secondary_genus_abundance

    # readlength
    Float? average_read_length = readlength.average_read_length
  }
}