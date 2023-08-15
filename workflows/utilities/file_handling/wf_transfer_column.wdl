version 1.0

import "../../../tasks/utilities/file_handling/task_file_handling.wdl" as file_handling
import "../../../tasks/task_versioning.wdl" as versioning

workflow transfer_column_content {
  input {
    Array[String] files_to_transfer
    String target_bucket
  }
  call file_handling.transfer_files {
    input:
      files_to_transfer = files_to_transfer,
      target_bucket = target_bucket
  }
  call versioning.version_capture {
    input:
  }
  output {
    String transfer_column_content_version = version_capture.phb_version
    String transfer_column_content_analysis_date = version_capture.date

    File transferred_files = transfer_files.transferred_files
  }
}