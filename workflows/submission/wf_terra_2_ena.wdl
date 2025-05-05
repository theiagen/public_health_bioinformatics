version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/data_handling/task_summarize_data.wdl" as task_summarize_data
import "../../tasks/utilities/submission/task_submit_ena_data.wdl" as task_submit_ena_data
import "../../tasks/utilities/submission/task_register_ena_samples.wdl" as task_register_ena_samples

workflow Terra_2_ENA {
  input {
    # Terra workspace and table information
    String terra_project_name
    String terra_workspace_name
    String terra_table_name
    Array[String] samples

    File? column_mappings # TSV file with column name mappings from Terra to ENA
    String sample_id_column # Column name containing sample IDs
    String study_accession # ENA study accession
    String ena_username # ENA username
    String ena_password # ENA password
    String sample_type # must be "prokaryotic_pathogen" or "virus_pathogen"
    Boolean allow_missing = false # Allow missing metadata for some samples
  }
  call task_summarize_data.summarize_data {
    input:
      sample_names = samples,
      terra_project = terra_project_name,
      terra_workspace = terra_workspace_name,
      terra_table = terra_table_name,
      output_prefix = "ena",
      id_column_name = sample_id_column
  }
  call task_register_ena_samples.register_ena_samples {
    input:
      metadata = summarize_data.filtered_metadata,
      ena_username = ena_username,
      ena_password = ena_password,
      study_accession = study_accession,
      column_mappings = column_mappings,
      sample_id_column = sample_id_column,
      sample_type = sample_type,
      allow_missing = allow_missing,
  }
  call task_submit_ena_data.submit_ena_data {
    input:
      ena_username = ena_username,
      ena_password = ena_password,
      metadata_accessions = register_ena_samples.metadata_accessions,
      study_accession = study_accession,
      column_mappings = column_mappings,
      sample_id_column = sample_id_column,
      allow_missing = allow_missing
  }
  call versioning.version_capture {
    input:
  }
  output {
    # Downloaded Terra table output
    File summarized_data = summarize_data.summarized_data
    File filtered_metadata = summarize_data.filtered_metadata
    # ENA registration output
    File ena_accessions = register_ena_samples.accessions
    File ena_metadata_accessions = register_ena_samples.metadata_accessions
    File ena_submission_summary = register_ena_samples.submission_summary
    File ena_submission_log = register_ena_samples.submission_log
    String ena_success = register_ena_samples.success
    # ENA data preparation output
    File prepped_ena_data = submit_ena_data.prepped_ena_data
    File ena_file_paths_json = submit_ena_data.file_paths_json
    File ena_excluded_samples = submit_ena_data.excluded_samples
    String ena_docker_image = submit_ena_data.docker_image
    # Versioning
    String Terra_2_ENA_version = version_capture.phb_version
    String Terra_2_ENA_analysis_date = version_capture.date
  }
}