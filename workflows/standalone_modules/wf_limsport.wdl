version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/data_export/task_download_terra_table.wdl" as download_table_task
import "../../tasks/utilities/data_handling/task_limsport.wdl" as limsport_task

workflow limsport_wf {
  meta {
    description: "Rename columns and QC rows of a TSV table to generate a LIMS-importable file using LIMSport"
  }
  input {
    File config
    Array[String] sample_names

    String terra_table_name
    String terra_workspace_name
    String terra_project_name
  }
  call download_table_task.download_terra_table {
    input:
      terra_table_name = terra_table_name,
      terra_workspace_name = terra_workspace_name,
      terra_project_name = terra_project_name
  }
  call limsport_task.limsport {
    input:
      input_table = download_terra_table.terra_table,
      config = config,
      sample_names = sample_names
  }
  call versioning.version_capture {
    input:
  }
  output {
    # PHB Version Captures
    String limsport_wf_version = version_capture.phb_version
    String limsport_wf_analysis_date = version_capture.date
    # LIMSport
    String limsport_version = limsport.limsport_version
    String limsport_docker = limsport.limsport_docker
    File limsport_table = limsport.limsport_table
    File limsport_qc_report = limsport.limsport_qc_report
  }
}
