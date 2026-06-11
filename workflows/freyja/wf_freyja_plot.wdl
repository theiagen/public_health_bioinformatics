version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/taxon_id/freyja/task_freyja_long_way.wdl" as freyja_long_format
import "../../tasks/taxon_id/freyja/task_freyja_microreact.wdl" as freyja_microreact
import "../../tasks/taxon_id/freyja/task_freyja_plot.wdl" as plot

workflow freyja_plot {
  input {
    Array[String] samplename
    Array[File] freyja_demixed
    Array[String]? freyja_lineages
    Array[String]? freyja_abundances
    Array[Float]? freyja_coverages
    Array[String]? collection_date
    Array[String]? collection_site
    Array[Float]? latitude
    Array[Float]? longitude
    String freyja_plot_name
    Int freyja_min_coverage = 60
    String freyja_long_format_docker = "us-docker.pkg.dev/general-theiagen/theiagen/freyja-microreact:1.0.2"
    String freyja_microreact_docker = "us-docker.pkg.dev/general-theiagen/theiagen/freyja-microreact:1.0.2"
  }
  String freyja_plot_name_updated = sub(freyja_plot_name, " ", "_")
  call plot.freyja_plot_task {
    input:
      samplename = samplename,
      freyja_demixed = freyja_demixed,
      collection_date = collection_date,
      freyja_plot_name = freyja_plot_name_updated,
      mincov = freyja_min_coverage
  }
  if (defined(freyja_lineages) && defined(freyja_abundances) && defined(freyja_coverages) && defined(collection_date) && defined(collection_site)) {
    call freyja_long_format.freyja_long_format_multi as freyja_long_format {
      input:
        samplenames = samplename,
        freyja_lineages = select_first([freyja_lineages]),
        freyja_abundances = select_first([freyja_abundances]),
        freyja_coverages = select_first([freyja_coverages]),
        collection_dates = collection_date,
        collection_sites = collection_site,
        latitudes = latitude,
        longitudes = longitude,
        mincov = freyja_min_coverage,
        docker = freyja_long_format_docker
    }
    call freyja_microreact.freyja_microreact {
      input:
        freyja_parsed_format_tsv = freyja_long_format.freyja_parsed_format_tsv,
        freyja_plot_name = freyja_plot_name_updated,
        docker = freyja_microreact_docker
    }
  }
  call versioning.version_capture {
    input:
  }
  output {
    # Freyja Plot Visualization
    File freyja_demixed_aggregate = freyja_plot_task.demixed_aggregate
    String? freyja_long_format_docker_used = freyja_long_format.freyja_long_format_docker
    String? freyja_microreact_docker_used = freyja_microreact.freyja_microreact_docker
    File? freyja_microreact_output = freyja_microreact.freyja_microreact_output
    File? freyja_parsed_format_tsv = freyja_long_format.freyja_parsed_format_tsv
    File freyja_plot = freyja_plot_task.freyja_plot
    File? freyja_plot_metadata = freyja_plot_task.freyja_plot_metadata
    String freyja_plot_version = freyja_plot_task.freyja_plot_version
    String freyja_plot_wf_analysis_date = version_capture.date
    String freyja_plot_wf_version = version_capture.phb_version
  }
}
