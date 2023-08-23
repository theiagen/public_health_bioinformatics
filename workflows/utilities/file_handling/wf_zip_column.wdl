version 1.0

import "../../../tasks/utilities/file_handling/task_file_handling.wdl" as file_handling
import "../../../tasks/task_versioning.wdl" as versioning

workflow zip_column_content {
  input {
    Array[File] files_to_zip
    String zipped_file_name
  }
  call file_handling.zip_files {
    input:
      files_to_zip = files_to_zip,
      zipped_file_name = zipped_file_name
	}
	call versioning.version_capture {
    input:
  }
  output {
    String zip_column_content_version = version_capture.phb_version
    String zip_column_content_analysis_date = version_capture.date

    File zipped_files = zip_files.zipped_files
  }
}