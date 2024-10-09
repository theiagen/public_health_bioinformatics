version 1.0

import "../../tasks/basecalling/task_dorado_basecall.wdl" as basecall_task

workflow dorado_basecalling_workflow {
  input {
    Array[File] input_files
    String dorado_model
  }

  call basecall_task.basecall {
    input:
      input_files = input_files,
      dorado_model = dorado_model
  }

  output {
    Array[File] basecalled_fastqs = basecall.basecalled_fastqs
    Array[File] logs = basecall.logs
  }
}
