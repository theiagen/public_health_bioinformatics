version 1.0

import "../../tasks/basecalling/task_dorado_basecall.wdl" as basecall_task

workflow dorado_basecalling_workflow {
  input {
    Array[File] input_files
    String dorado_model
    String kit_name
  }

  call basecall_task.basecall {
    input:
      input_files = input_files,
      dorado_model = dorado_model,
      kit_name = kit_name
  }

  output {
    Array[File] combined_fastqs = basecall.combined_fastqs
    Array[File] logs = basecall.logs
  }
}
