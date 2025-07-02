version 1.0

import "../../tasks/species_typing/escherichia_shigella/task_shigapass.wdl" as shigapass
import "../../tasks/task_versioning.wdl" as versioning

workflow shigapass_workflow {
  input {
    Array[File] assembly_fasta
    Array[String] samplename
  }
  call shigapass.shigapass_many as shigapass_task {
    input:
      assembly_fasta = assembly_fasta,
      samplename = samplename
  }
  call versioning.version_capture {
    input:
  }
  output {
    # Version Capture
    String shigapass_wf_version = version_capture.phb_version
    String shigapass_wf_analysis_date = version_capture.date
    # shigapass outputs
    # NOTE: may end up removing the shigapass_summary and shigapass_flexneri_summary outputs since they are SEMICOLON delmited sheets
    # and not TSV. Keeping for now to ensure the conversion from semi-colon delimited to TSV is working as expected.
    File shigapass_summary = shigapass_task.shigapass_summary
    File shigapass_summary_tsv = shigapass_task.shigapass_summary_tsv
    File? shigapass_flexneri_summary = shigapass_task.shigapass_flexneri_summary
    File? shigapass_flexneri_summary_tsv = shigapass_task.shigapass_flexneri_summary_tsv
    String shigapass_version = shigapass_task.shigapass_version
    String shigapass_docker = shigapass_task.shigapass_docker
  }
}