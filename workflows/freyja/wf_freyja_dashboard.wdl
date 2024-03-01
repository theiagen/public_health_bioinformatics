version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/taxon_id/freyja/task_freyja_dashboard.wdl" as freyja_dash

workflow freyja_dashboard {
  input {
    Array[String] samplename
    Array[File] freyja_demixed
    Array[String] collection_date
    Array[String] viral_load
    String freyja_dashboard_title
    File? dashboard_intro_text
  }
  call freyja_dash.freyja_dashboard_task {
    input:
      samplename = samplename,
      freyja_demixed = freyja_demixed,
      collection_date = collection_date,
      viral_load = viral_load,
      freyja_dashboard_title = freyja_dashboard_title,
      dashboard_intro_text = dashboard_intro_text,
  }
  call versioning.version_capture {
    input:
  }
  output {
    # Version Capture
    String freyja_dashboard_wf_version = version_capture.phb_version
    String freyja_dashboard_wf_analysis_date = version_capture.date
    # Freyja Dashboard Visualization
    String freyja_dashboard_version = freyja_dashboard_task.freyja_dashboard_version
    File freyja_dasbhoard = freyja_dashboard_task.freyja_dasbhoard
    File freyja_demixed_aggregate = freyja_dashboard_task.freyja_demixed_aggregate
    File freyja_dashboard_metadata = freyja_dashboard_task.freyja_dashboard_metadata
  }
}