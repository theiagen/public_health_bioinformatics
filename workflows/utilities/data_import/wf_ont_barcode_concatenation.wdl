version 1.0

import "../../../tasks/utilities/file_handling/task_cat_ont_barcodes.wdl" as cat_ont_barcodes_task
import "../../../tasks/utilities/data_import/task_create_terra_table.wdl" as create_terra_table_task

workflow ont_barcode_concatenation {
  input {
    String input_bucket_path
    String output_bucket_path

    String terra_table_name
    String terra_project
    String terra_workspace

    File? barcode_renaming_file
    String file_extension = ".fastq.gz"
  }
  call cat_ont_barcodes_task.cat_ont_barcodes {
    input:
      input_bucket_path = input_bucket_path,
      output_bucket_path = output_bucket_path,
      file_extension = file_extension,
      barcode_renaming_file = barcode_renaming_file
  } 
  if (defined(cat_ont_barcodes.concatenation_log)) {
    call create_terra_table_task.create_terra_table {
      input:
        new_table_name = terra_table_name,
        terra_project = terra_project,
        terra_workspace = terra_workspace,
        data_location_path = output_bucket_path,
        paired_end = false,
        assembly_data = false,
        file_ending = file_extension,
        responsible_workflow = "ONT_Barcode_Concatenation_PHB"
    }
  }
  output {
  }
}