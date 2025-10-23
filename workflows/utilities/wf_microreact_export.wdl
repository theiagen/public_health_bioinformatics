version 1.0

import "../../tasks/utilities/data_export/task_download_terra_table.wdl" as task_download_terra_table
import "../../tasks/utilities/data_export/task_microreact_export.wdl" as microreact_export

workflow wf_microreact_export {
  input {
    String project_name
    String id_column
    String terra_table_name
    String terra_workspace_name
    String terra_project_name
    String? project_url
    String? access_token
    File? metadata_file
    Array[File]? tree_files
    Array[String]? metadata_columns
    Boolean update_project = false
    Boolean restricted_access = true
    Boolean download_terra_table = true
  }
  if (download_terra_table) {
    call task_download_terra_table.download_terra_table as download {
      input:
        terra_table_name = terra_table_name,
        terra_workspace_name = terra_workspace_name,
        terra_project_name = terra_project_name
    }
  }
  call microreact_export.microreact_export as microreact{
    input:
      project_name = project_name,
      metadata_tsv = select_first([download.terra_table, metadata_file]),
      tree_files = tree_files,
      metadata_columns = metadata_columns,
      update_project = update_project,
      project_url = project_url,
      restricted_access = restricted_access,
      access_token = access_token,
      id_column = id_column
  }
  output {
    File microreact_json = microreact.microreact_json
    File? microreact_api_response = microreact.microreact_api_response
  }
}