version 1.0

import "../../../tasks/task_versioning.wdl" as versioning
import "../../../tasks/utilities/file_handling/task_cat_files.wdl" as file_handling

workflow concatenate_column_content {
  input {
    Array[File] files_to_cat
    String concatenated_file_name
    String concatenated_file_name_updated = sub(concatenated_file_name, " ", "_")
  }
  call file_handling.cat_files {
    input:
      files_to_cat = files_to_cat,
      concatenated_file_name = concatenated_file_name_updated
  }
  call versioning.version_capture {
    input:
  }
  output {
    String concatenate_column_content_version = version_capture.phb_version
    String concatenate_column_content_analysis_date = version_capture.date

    File concatenated_files = cat_files.concatenated_files
  }
}