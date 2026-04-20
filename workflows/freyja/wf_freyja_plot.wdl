version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/taxon_id/freyja/task_freyja_plot.wdl" as plot
import "../../tasks/taxon_id/freyja/task_freyja_long_way.wdl" as freyja_long_way
import "../../tasks/taxon_id/freyja/task_freyja_microreact.wdl" as freyja_microreact

workflow freyja_plot {
  input {
    Array[String] samplename
    Array[File] freyja_demixed
    Array[String]? freyja_lineages
    Array[String]? freyja_abundances
    Array[String]? collection_date
    Array[String]? collection_site
    Array[Float]? latitude
    Array[Float]? longitude
    String freyja_plot_name
  }
  String freyja_plot_name_updated = sub(freyja_plot_name, " ", "_")
  call plot.freyja_plot_task {
    input:
      samplename = samplename,
      freyja_demixed = freyja_demixed,
      collection_date = collection_date,
      freyja_plot_name = freyja_plot_name_updated
  }
  if (defined(freyja_lineages) && defined(freyja_abundances) && defined(collection_date) && defined(collection_site)) {
    call freyja_long_way.freyja_long_way_multi {
      input:
        samplenames       = samplename,
        freyja_lineages   = select_first([freyja_lineages]),
        freyja_abundances = select_first([freyja_abundances]),
        collection_dates  = collection_date,
        collection_sites  = collection_site,
        latitudes         = latitude,
        longitudes        = longitude
    }
    call freyja_microreact.freyja_microreact {
      input:
        freyja_long_format_tsv = freyja_long_way_multi.freyja_long_format,
        freyja_plot_name = freyja_plot_name_updated
    }
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
    # Freyja Long Way
    File? freyja_long_format = freyja_long_way_multi.freyja_long_format
    # Freyja Microreact
    File? freyja_microreact_output = freyja_microreact.freyja_microreact_output
  }
}