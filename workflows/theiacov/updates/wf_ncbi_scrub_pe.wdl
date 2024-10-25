version 1.0

import "../../../tasks/quality_control/read_filtering/task_ncbi_scrub.wdl" as ncbi_scrub
import "../../../tasks/task_versioning.wdl" as versioning
import "../../../tasks/taxon_id/contamination/task_kraken2.wdl" as kraken

workflow dehost_pe {
  input {
    String samplename
    File read1
    File read2
    File kraken2_db = "gs://theiagen-large-public-files-rp/terra/databases/kraken2/kraken2_humanGRCh38_viralRefSeq_20240828.tar.gz"

  }
  call ncbi_scrub.ncbi_scrub_pe {
    input:
      samplename = samplename,
      read1 = read1,
      read2 = read2
  }
  call kraken.kraken2_standalone as kraken2 {
    input:
      samplename = samplename,
      read1 = ncbi_scrub_pe.read1_dehosted,
      read2 = ncbi_scrub_pe.read2_dehosted,
      kraken2_db = kraken2_db,
      target_organism = "Severe acute respiratory syndrome coronavirus 2"
  }
  call versioning.version_capture {
    input:
  }
  output {
    String ncbi_scrub_pe_version = version_capture.phb_version
    String ncbi_scrub_pe_analysis_date = version_capture.date
    File read1_dehosted = ncbi_scrub_pe.read1_dehosted
    File read2_dehosted = ncbi_scrub_pe.read2_dehosted
    Int ncbi_scrub_human_spots_removed = ncbi_scrub_pe.human_spots_removed
    String ncbi_scrub_docker = ncbi_scrub_pe.ncbi_scrub_docker
    Float kraken2_human_dehosted = kraken2.kraken2_percent_human
    Float kraken2_sc2_dehosted = kraken2.kraken2_percent_sc2
    File kraken2_report_dehosted = kraken2.kraken2_report
    String kraken2_version_dehosted = kraken2.kraken2_version
  }
}