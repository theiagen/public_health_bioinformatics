version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/taxon_id/freyja/task_freyja_plot.wdl" as plot

workflow freyja_plot {
  input {
    Array[String] samplename
    Array[File] freyja_demixed
    Array[String]? collection_date
    String freyja_plot_name
  }
  call plot.freyja_plot_task {
    input:
      samplename = samplename,
      freyja_demixed = freyja_demixed,
      collection_date = collection_date,
      freyja_plot_name = freyja_plot_name
  }
  call versioning.version_capture {
    input:
  }
  output {
    # Version Capture
    String freyja_plot_wf_version = version_capture.phb_version
    String freyja_plot_wf_analysis_date = version_capture.date
    # Freyja Plot Visualization
    String freyja_plot_version = freyja_plot_task.freyja_plot_version
    File freyja_plot = freyja_plot_task.freyja_plot
    File freyja_demixed_aggregate = freyja_plot_task.demixed_aggregate
    File? freyja_plot_metadata = freyja_plot_task.freyja_plot_metadata
  }
}