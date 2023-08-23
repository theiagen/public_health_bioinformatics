version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/data_export/task_export_two_tsvs.wdl" as export
import "../../tasks/utilities/data_handling/task_validate.wdl" as validate_task

workflow theiavalidate {
  input {
    String table1 # consider name change
    String table2
    String columns_to_compare # comma separated
    String output_prefix

    # table 1 location
    String terra_workspace1_name
    String terra_project1_name

    # table 2 location (optional)
    String? terra_workspace2_name
    String? terra_project2_name

    # criteria to test (optional)
    File? validation_criteria_tsv
  }
  call export.export_two_tsvs {
    input:
      datatable1 = table1,
      terra_workspace1 = terra_workspace1_name,
      terra_project1 = terra_project1_name,
      datatable2 = table2,
      terra_workspace2 = terra_workspace2_name,
      terra_project2 = terra_project2_name
  }
  if (!export_two_tsvs.same_table_length) {
    String validation_failure = "Input tables were not of same length; validation not performed"
  }
  if (export_two_tsvs.same_table_length) {
    String validation_attempted = "Validation attempted"
    call validate_task.validate{
      input:
        datatable1 = table1,
        datatable1_tsv = export_two_tsvs.datatable1_tsv,
        datatable2 = table2,
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
    File? validation_report = validate.pdf_report
    File? validation_differences_table = validate.excel_report
    File? input_table1 = validate.input_table1
    File? input_table2 = validate.input_table2
    String validation_status = select_first([validation_failure, validation_attempted])
  }
}
