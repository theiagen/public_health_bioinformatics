version 1.0

import "../../tasks/basecalling/task_dorado_basecall.wdl" as basecall_task
import "../../tasks/basecalling/task_samtools_convert.wdl" as samtools_convert_task
import "../../tasks/basecalling/task_dorado_demux.wdl" as dorado_demux_task

workflow dorado_basecalling_workflow {
  input {
    Array[File] input_files
    String dorado_model
    String kit_name
  }

  call basecall_task.basecall as basecall_step {
    input:
      input_files = input_files,
      dorado_model = dorado_model,
      kit_name = kit_name,
  }

  call samtools_convert_task.samtools_convert as samtools_convert_step {
    input:
      sam_files = basecall_step.sam_files
  }

  call dorado_demux_task.dorado_demux as dorado_demux_step {
    input:
      bam_files = samtools_convert_step.bam_files,
      kit_name = kit_name
  }

  output {
    Array[File] fastq_files = dorado_demux_step.fastq_files
  }
}
