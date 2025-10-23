version 1.0

import "../../tasks/utilities/data_import/task_lod_table_prep.wdl" as lod_table_prep_task

workflow lod_prep {
  input {
    File lod_config_yaml

  }
  call lod_table_prep_task.lod_table_prep as lod_table_prep {
    input:
      lod_config_yaml = lod_config_yaml
  }
  output {
    File lod_table_tsv = lod_table_prep.lod_table_tsv
  } 
}