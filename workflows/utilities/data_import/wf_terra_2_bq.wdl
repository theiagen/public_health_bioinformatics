version 1.0

import "../../../tasks/utilities/data_import/task_terra_2_bq.wdl" as terra_2_bq_task

workflow terra_2_bq {
    input {
      Array[String]  terra_projects
      Array[String]  workspace_names
      Array[String]  table_names # name of Terra data table, used to download data table via firecloud API
      Array[String]  table_ids # ID column or the 1st column is not really necessary - we are not writing back to the data table. Removing to reduce variable name confusion
      Array[String]  gcs_uris
      Array[String]  output_filename_prefixs
    }
    call terra_2_bq_task.terra_to_bigquery {
      input:
        terra_projects=terra_projects,
        workspace_names=workspace_names,
        table_names=table_names,
        table_ids=table_ids,
        gcs_uri_prefixs=gcs_uris,
        output_filename_prefix=output_filename_prefixs
    }
}

