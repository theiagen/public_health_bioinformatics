version 1.0

import "../../../tasks/quality_control/advanced_metrics/task_vadr.wdl" as vadr_task
import "../../../tasks/quality_control/basic_statistics/task_consensus_qc.wdl" as consensus_qc_task
import "../../utilities/wf_organism_parameters.wdl" as set_organism_defaults
import "../../../tasks/task_versioning.wdl" as versioning

workflow vadr_update {
  input {
    File genome_fasta
    String organism = "sars-cov-2" # options: "sars-cov-2" "MPXV" "WNV" "flu" "rsv_a" "rsv_b" "measles" "mumps" "rubella"
    # vadr parameters
    Int? vadr_max_length
    Int? vadr_skip_length
    String? vadr_opts
    File? vadr_model_file
    Int? vadr_memory
  }
  call set_organism_defaults.organism_parameters {
    input:
      organism = organism,
      vadr_max_length = vadr_max_length,
      vadr_skip_length = vadr_skip_length,
      vadr_options = vadr_opts,
      vadr_model = vadr_model_file,
      vadr_mem = vadr_memory
  }
  call consensus_qc_task.consensus_qc {
    input:
      assembly_fasta = genome_fasta
  }
  if (organism_parameters.standardized_organism == "sars-cov-2" || organism_parameters.standardized_organism == "MPXV" || organism_parameters.standardized_organism == "rsv_a" || organism_parameters.standardized_organism == "rsv_b" || organism_parameters.standardized_organism == "WNV" || organism_parameters.standardized_organism == "flu" || organism_parameters.standardized_organism == "mumps" || organism_parameters.standardized_organism == "rubella" || organism_parameters.standardized_organism == "measles") {
    call vadr_task.vadr {
      input:
        genome_fasta = genome_fasta,
        assembly_length_unambiguous = consensus_qc.number_ATCG,
        max_length = organism_parameters.vadr_maxlength,
        vadr_opts = organism_parameters.vadr_opts,
        vadr_model_file = organism_parameters.vadr_model_file,
        skip_length = organism_parameters.vadr_skiplength,
        memory = organism_parameters.vadr_memory
    }
  }
  call versioning.version_capture {
    input:
  }
  output {
    # Version Capture
    String vadr_update_version = version_capture.phb_version
    String vadr_update_analysis_date = version_capture.date
    # VADR Annotation QC
    File? vadr_alerts_list = vadr.alerts_list
    File? vadr_feature_tbl_pass = vadr.feature_tbl_pass
    File? vadr_feature_tbl_fail = vadr.feature_tbl_fail
    File? vadr_classification_summary_file = vadr.classification_summary_file
    File? vadr_all_outputs_tar_gz = vadr.outputs_tgz
    String? vadr_num_alerts = vadr.num_alerts
    String? vadr_docker = vadr.vadr_docker
    File? vadr_fastas_zip_archive = vadr.vadr_fastas_zip_archive
  }
}