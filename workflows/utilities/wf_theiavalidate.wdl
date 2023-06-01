version 1.0

import "../../tasks/utilities/task_validate.wdl" as validate
import "../../tasks/task_versioning.wdl" as versioning

workflow theiavalidate {
  input {
    String old_version_table # consider name change
    String new_version_table
    String terra_workspace_name
    String terra_project_name
    File? validation_criteria_tsv
    String columns_to_compare
    String output_prefix
  }
  call validate.export_two_tsvs {
    input:
      terra_project = terra_project_name,
      terra_workspace = terra_workspace_name,
      datatable1 = old_version_table,
      datatable2 = new_version_table
  }
  if (!export_two_tsvs.same_table_length) {
    String validation_failure = "Input tables were not of same length; validation not performed"
  }
  if (export_two_tsvs.same_table_length) {
    String validation_attempted = "Validation attempted"
    call validate.compare_two_tsvs{
      input:
        datatable1 = old_version_table,
        datatable1_tsv = export_two_tsvs.datatable1_tsv,
        datatable2 = new_version_table,
        datatable2_tsv = export_two_tsvs.datatable2_tsv,
        validation_criteria_tsv = validation_criteria_tsv,
        columns_to_compare = columns_to_compare,
        output_prefix = output_prefix
    }
  }
  call versioning.version_capture {
    input:
  }
  output {
    String theiavalidate_version = version_capture.phb_version
    String theiavalidate_date = version_capture.date
    File? validation_report = compare_two_tsvs.pdf_report
    File? validation_differences_table = compare_two_tsvs.excel_report
    String validation_status = select_first([validation_failure, validation_attempted])
  }
}
