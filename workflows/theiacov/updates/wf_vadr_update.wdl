version 1.0

import "../../../tasks/quality_control/task_vadr.wdl" as vadr_task
import "../../../tasks/task_versioning.wdl" as versioning

workflow vadr_update {
  input {
    File genome_fasta
    String docker
    Int assembly_length_unambiguous
  }
  call vadr_task.vadr {
    input:
      genome_fasta = genome_fasta,
      docker = docker,
      assembly_length_unambiguous = assembly_length_unambiguous
  }
  call versioning.version_capture{
    input:
  }
  output {
    # Version Capture
    String vadr_update_version = version_capture.phb_version
    String vadr_update_analysis_date = version_capture.date
    # VADR Annotation QC
    File? vadr_alerts_list = vadr.alerts_list
    String vadr_num_alerts = vadr.num_alerts
    String vadr_docker = vadr.vadr_docker
    File? vadr_fastas_zip_archive = vadr.vadr_fastas_zip_archive
  }
}