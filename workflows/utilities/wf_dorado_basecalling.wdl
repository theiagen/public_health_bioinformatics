version 1.0

import "../../tasks/basecalling/dorado_basecall.wdl" as dorado_basecall

workflow dorado_basecalling_workflow {
  input {
    Array[File] input_files
    Array[String] sample_names
    String dorado_model
    String output_prefix
  }

  call dorado_basecall {
    input:
      input_files = input_files,
      sample_names = sample_names,
      dorado_model = dorado_model,
      output_prefix = output_prefix
  }

  output {
    Array[File] basecalled_fastqs = dorado_basecall.basecalled_fastqs
    Array[File] logs = dorado_basecall.logs
  }
}
