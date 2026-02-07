version 1.0

import "../../../tasks/quality_control/read_filtering/task_ncbi_scrub.wdl" as ncbi_scrub
import "../../../tasks/task_versioning.wdl" as versioning
import "../../../tasks/taxon_id/contamination/task_kraken2.wdl" as kraken

workflow dehost_se {
  input {
    String samplename
    File read1
    String target_organism = "Severe acute respiratory syndrome coronavirus 2"
  }
  call ncbi_scrub.ncbi_scrub_se {
    input:
      samplename = samplename,
      read1 = read1
  }
  call kraken.kraken2_theiacov as kraken2 {
    input:
      samplename = samplename,
      read1 = ncbi_scrub_se.read1_dehosted,
      target_organism = target_organism
  }
  call versioning.version_capture {
    input:
  }
  output {
    String ncbi_scrub_se_version = version_capture.phb_version
    String ncbi_scrub_se_analysis_date = version_capture.date
    File read1_dehosted = ncbi_scrub_se.read1_dehosted
    String ncbi_scrub_docker = ncbi_scrub_se.ncbi_scrub_docker
    Int ncbi_scrub_human_spots_removed = ncbi_scrub_se.human_spots_removed
    Float kraken_human_dehosted = kraken2.percent_human
    String kraken_sc2_dehosted = kraken2.percent_sc2
    String kraken_version_dehosted = kraken2.kraken2_version
    File kraken_report_dehosted = kraken2.kraken_report
  }
}