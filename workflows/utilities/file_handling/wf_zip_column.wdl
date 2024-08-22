version 1.0

import "../../../tasks/task_versioning.wdl" as versioning
import "../../../tasks/utilities/file_handling/task_zip_files.wdl" as file_handling

workflow zip_column_content {
  input {
    Array[File] files_to_zip
    String zipped_file_name
  }
  String zipped_file_name_updated = sub(zipped_file_name, " ", "_")
  call file_handling.zip_files {
    input:
      files_to_zip = files_to_zip,
      zipped_file_name = zipped_file_name_updated
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