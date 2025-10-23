version 1.0

import "../tasks/utilities/task_download_terra_table.wdl" as task_download_terra_table
import "../tasks/utilities/data_export/task_microreact_export.wdl" as microreact_export

workflow wf_microreact_export {
  input {
    String terra_table_name
    String terra_workspace_name
    String terra_project_name
    Array[String] tree_files
    Array[String] metadata_columns
  }
  call task_download_terra_table.download_terra_table {
    input:
      terra_table_name = terra_table_name,
      terra_workspace_name = terra_workspace_name,
      terra_project_name = terra_project_name
  }
  call microreact_export.microreact_export {
    input:
      metadata_tsv = terrain_download_terra_table.terra_table,
      tree_files = tree_files,
      metadata_columns = metadata_columns,
  }
  output {

  }
}