version 1.0

import "../../tasks/basecalling/task_dorado_basecall.wdl" as basecall_task

workflow dorado_basecalling_workflow {
  input {
    Array[File] input_files
    Array[String] sample_names
    String dorado_model
    String output_prefix
  }

  call basecall_task.dorado_basecall {  # Use the new alias here
    input:
      input_files = input_files,
      sample_names = sample_names,
      dorado_model = dorado_model,
      output_prefix = output_prefix
  }

  output {
    Array[File] basecalled_fastqs = basecall_task.dorado_basecall.basecalled_fastqs  # Use the new alias here
    Array[File] logs = basecall_task.dorado_basecall.logs  # Use the new alias here
  }
}
