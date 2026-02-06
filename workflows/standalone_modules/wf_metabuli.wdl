version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/taxon_id/contamination/task_metabuli.wdl" as metabuli
import "../../tasks/taxon_id/contamination/task_krona.wdl" as krona_task

workflow metabuli_wf {
  meta {
    description: "Classify paired-end reads using Metabuli"
  }
  input {
    String samplename
    File read1
    File? read2
    File metabuli_db
  }
  call metabuli.metabuli as metabuli {
    input:
      samplename = samplename,
      read1 = read1,
      read2? = read2,
      metabuli_db = metabuli_db
  }
  call krona_task.krona as krona {
    input:
      metabuli_report = metabuli.metabuli_report,
      samplename = samplename
  }
  call versioning.version_capture {
    input:
  }
  output {
    # PHB Version Captures
    String metabuli_wf_version = version_capture.phb_version
    String metabuli_wf_analysis_date = version_capture.date
    # Metabuli
    String metabuli_version = metabuli.metabuli_version
    String metabuli_docker = metabuli.metabuli_docker
    File metabuli_report = metabuli.metabuli_report
    File metabuli_classified_report = metabuli.metabuli_classified
    File metabuli_classified_read1 = metabuli.metabuli_read1_extract
    # Krona outputs
    String krona_version = krona.krona_version
    String krona_docker = krona.krona_docker
    File krona_html = krona.krona_html
  }
}
