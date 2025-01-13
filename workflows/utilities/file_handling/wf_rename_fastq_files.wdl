version 1.0

import "../../../tasks/task_versioning.wdl" as versioning
import "../../../tasks/utilities/file_handling/task_rename_files.wdl" as rename_files_task

workflow rename_fastq_files {
  input {
    File read1
    File? read2
    String new_filename
  }
  String new_filename_updated = sub(new_filename, " ", "_")
  if (defined(read2)) {
    call rename_files_task.rename_PE_files {
      input:
        read1 = read1,
        read2 = select_first([read2]),
        new_filename = new_filename_updated
    }
  } 
  if (!defined(read2)) {
    call rename_files_task.rename_SE_files {
      input:
        read1 = read1,
        new_filename = new_filename_updated
    }
  }
  call versioning.version_capture {
    input:
  }
  output {
    String rename_fastq_files_version = version_capture.phb_version
    String rename_fastq_files_analysis_date = version_capture.date
    File read1_renamed = select_first([rename_PE_files.read1_renamed, rename_SE_files.read1_renamed])
    File? read2_renamed = rename_PE_files.read2_renamed
  }
}