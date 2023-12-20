version 1.0

import "../../tasks/taxon_id/task_kraken2.wdl" as kraken2
import "../../tasks/taxon_id/task_krona.wdl" as krona_task
import "../../tasks/task_versioning.wdl" as versioning

workflow kraken2_se_wf {
  meta {
    description: "Classify single-end reads using Kraken2"
  }
  input {
    String  samplename
    File    read1
    File    kraken2_db
  }
  call kraken2.kraken2_standalone as kraken2_se {
    input:
      samplename = samplename,
      read1 = read1,
      kraken2_db = kraken2_db
  }
  call krona_task.krona as krona {
    input:
      kraken2_report = kraken2_se.kraken2_report,
      samplename = samplename
  }
  call versioning.version_capture{
    input:
  }
  output {
    # PHBG Version Captures
    String kraken2_se_wf_version = version_capture.phb_version
    String kraken2_se_wf_analysis_date = version_capture.date
    # Kraken2
    String kraken2_version = kraken2_se.kraken2_version
    String kraken2_docker = kraken2_se.kraken2_docker
    File kraken2_report = kraken2_se.kraken2_report
    File kraken2_classified_report = kraken2_se.kraken2_classified_report
    File kraken2_unclassified_read1 = kraken2_se.kraken2_unclassified_read1
    File kraken2_classified_read1 = kraken2_se.kraken2_classified_read1
    # Krona outputs
    String krona_version = krona.krona_version
    String krona_docker = krona.krona_docker
    File krona_html = krona.krona_html
  }
}
