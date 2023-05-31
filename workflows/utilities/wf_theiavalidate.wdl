version 1.0

import "../../tasks/utilities/task_validate.wdl" as validate
import "../../tasks/task_versioning.wdl" as versioning

workflow theiavalidate {
  input {
    String old_version_table # consider name change
    String new_version_table
    String terra_workspace_name
    String terra_project_name
    File validation_criteria_tsv
    String columns_to_compare
  }
  call validate.export_two_tsvs {
    input:
      terra_project = terra_project_name,
      terra_workspace = terra_workspace_name,
      datatable1 = old_version_table,
      datatable2 = new_version_table
  }
  call versioning.version_capture {

  }
  output {

  }
}
