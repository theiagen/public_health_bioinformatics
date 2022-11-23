version 1.0

import "../tasks/utilities/task_rasusa.wdl" as rasusa
import "../tasks/task_versioning.wdl" as versioning

workflow rasusa_workflow {
  input {
    File read1
    File? read2
    String samplename
    Float coverage
    String genome_size
  }
  call rasusa.rasusa as rasusa_task {
    input:
      read1 = read1,
      read2 = read2,
      samplename = samplename,
      genome_size = genome_size,
      coverage = coverage
  }
  call versioning.version_capture{
    input:
  }
  output {
    String rasusa_wf_version = version_capture.phbg_version
    String rasusa_wf_analysis_date = version_capture.date

    String rasusa_version = rasusa_task.rasusa_version
    File read1_subsampled = rasusa_task.read1_subsampled
    File? read2_subsampled = rasusa_task.read2_subsampled
  }
}