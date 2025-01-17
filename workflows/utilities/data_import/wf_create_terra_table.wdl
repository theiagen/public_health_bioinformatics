version 1.0

import "../../../tasks/utilities/data_import/task_create_terra_table.wdl" as make_table_task

workflow create_terra_table {
  input {
    String new_table_name
    String data_location_path
    Boolean paired_end = false
    Boolean assembly_data = false

    String file_ending = ".fastq.gz"

    String terra_project
    String terra_workspace
  }
  call make_table_task.create_terra_table as make_table {
    input:
      new_table_name = new_table_name,
      data_location_path = data_location_path,
      paired_end = paired_end,
      assembly_data = assembly_data,
      file_ending = file_ending,
      terra_project = terra_project,
      terra_workspace = terra_workspace
  }
  output {
  }
}