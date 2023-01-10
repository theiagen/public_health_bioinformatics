version 1.0

import "../tasks/phylogenetic_inference/task_lyveset.wdl" as lyvset
import "../tasks/task_versioning.wdl" as versioning

workflow lyvset_workflow {
  input {
    Array[File] read1
    Array[File] read2
    String samplename
    File reference_genome
  }
  call lyveset.lyvset as lyvset_task {
    input:
      read1 = read1,
      read2 = read2,
      samplename = samplename,
      reference_genome = reference_genome
  }
  call versioning.version_capture{
    input:
  }
  output {
    String lyveset_wf_version = version_capture.phb_version
    String lyvset_wf_analysis_date = version_capture.date

    String rasusa_version = rasusa_task.rasusa_version
    File read1_subsampled = rasusa_task.read1_subsampled
    File? read2_subsampled = rasusa_task.read2_subsampled
  }
}