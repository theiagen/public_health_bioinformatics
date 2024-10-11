version 1.0

import "../../tasks/assembly/task_flye.wdl" as flye_task

workflow flye_consensus {
  meta {
    description: "This workflow runs either flye, miniasm, or raven to generate a consensus genome assembly from long reads."
  }
  input {
    File read1
    String samplename
  }
  call flye_task.flye {
    input:
      read1 = read1,
      samplename = samplename
  }
  # placeholder
}