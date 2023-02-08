version 1.0

import "../tasks/task_read_clean.wdl" as read_clean
import "../tasks/task_taxonID.wdl" as taxonID
import "../tasks/task_versioning.wdl" as versioning

workflow dehost_se {
  input {
    String samplename
    File reads
  }
  call read_clean.ncbi_scrub_se {
    input:
      samplename = samplename,
      read1 = reads
  }
  call taxonID.kraken2 {
    input:
      samplename = samplename,
      read1 = ncbi_scrub_se.read1_dehosted
  }
  call versioning.version_capture {
    input:
  }
  output {
    String ncbi_scrub_se_version = version_capture.phvg_version
    String ncbi_scrub_se_analysis_date = version_capture.date
    File reads_dehosted = ncbi_scrub_se.read1_dehosted
    String ncbi_scrub_docker = ncbi_scrub_se.ncbi_scrub_docker
    Int human_spots_removed = ncbi_scrub_se.read1_human_spots_removed
    Float kraken_human_dehosted = kraken2.percent_human
    Float kraken_sc2_dehosted = kraken2.percent_sc2
    String kraken_version_dehosted = kraken2.version
    File kraken_report_dehosted = kraken2.kraken_report
  }
}