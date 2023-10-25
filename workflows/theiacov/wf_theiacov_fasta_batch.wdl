version 1.0

import "../../tasks/species_typing/task_pangolin.wdl" as pangolin_task
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade_task
import "../../tasks/utilities/task_file_handling.wdl" as concatenate

workflow theiacov_fasta_batch {
  meta {
    description: "TheiaCoV_FASTA for multiple samples"
  }
  input {
    Array[String] samplenames
    Array[File] assembly_fastas
    String organism = "sars-cov-2"
    # sequencing values
    String seq_method
    String input_assembly_method
    # nextclade inputs
    String nextclade_dataset_reference = "MN908947"
    String nextclade_dataset_tag = "2023-09-21T12:00:00Z"
    String? nextclade_dataset_name
  }
  call concatenate.cat_files {
    input: 
      files_to_cat = assembly_fastas,
      concatenated_file_name = "concatenated_assemblies.fasta"
  }
  call 
}