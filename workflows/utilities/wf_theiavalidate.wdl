version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/data_handling/task_validate.wdl" as validate

workflow theiavalidate {
  input {
    String table1_name # consider name change
    String table2_name
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

    # column translation (optional)
    File? column_translation_tsv
  }
  call validate.export_two_tsvs {
    input:
      datatable1 = table1_name,
      terra_workspace1 = terra_workspace1_name,
      terra_project1 = terra_project1_name,
      datatable2 = table2_name,
      terra_workspace2 = terra_workspace2_name,
      terra_project2 = terra_project2_name
  }
  if (!export_two_tsvs.same_table_length) {
    String validation_failure = "Input tables were not of same length; validation not performed"
  }
  if (export_two_tsvs.same_table_length) {
    String validation_attempted = "Validation attempted"
    call validate.theiavalidate as compare_two_tsvs {
      input:
        datatable1_tsv = export_two_tsvs.datatable1_tsv,
        datatable2_tsv = export_two_tsvs.datatable2_tsv,
        columns_to_compare = columns_to_compare,
        output_prefix = output_prefix,
        validation_criteria_tsv = validation_criteria_tsv,
        column_translation_tsv = column_translation_tsv
    }
  }
  call versioning.version_capture {
    input:
  }
  output {
    String theiavalidate_wf_version = version_capture.phb_version
    String theiavalidate_date = version_capture.date
    String? theiavalidate_version = compare_two_tsvs.theiavalidate_version
    File? theiavalidate_report = compare_two_tsvs.summary_pdf_report
    File? theiavalidate_exact_differences = compare_two_tsvs.exact_differences
    File? theiavalidate_criteria_differences = compare_two_tsvs.validation_criteria_differences
    File? theiavalidate_filtered_input_table1 = compare_two_tsvs.filtered_input_table1
    File? theiavalidate_filtered_input_table2 = compare_two_tsvs.filtered_input_table2
    String theiavalidate_status = select_first([validation_failure, validation_attempted])
    Array[File]? theiavalidate_diffs = compare_two_tsvs.diffs
  }
}
