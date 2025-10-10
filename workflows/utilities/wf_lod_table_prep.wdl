version 1.0

import "../../tasks/utilities/data_import/task_lod_table_prep.wdl" as lod_table_prep_task

workflow lod_prep {
  input {
    String input_table_name
    String workspace_name
    String project_name

    String output_table_name

    String read1_column_name = "read1"
    String read2_column_name = "read2"
    String taxon_column_name = "gambit_predicted_taxon"

    Array[Int] downsampling_levels
    Array[String] expected_genes
    Array[String] expected_alleles

  }
  call lod_table_prep_task.lod_table_prep as lod_table_prep {
    input:
      input_table_name = input_table_name,
      workspace_name = workspace_name,
      project_name = project_name,
      output_table_name = output_table_name,
      read1_column_name = read1_column_name,
      read2_column_name = read2_column_name,
      taxon_column_name = taxon_column_name,
      downsampling_levels = downsampling_levels,
      expected_genes = expected_genes,
      expected_alleles = expected_alleles
  }
  output {
    File lod_table_tsv = lod_table_prep.lod_table_tsv
  } 
}