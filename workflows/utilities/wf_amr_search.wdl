version 1.0

import "../../tasks/gene_typing/drug_resistance/task_amr_search.wdl" as run_amr_search

workflow amr_search_workflow {
  input {
    File input_fasta
    String amr_search_database
    String samplename
  }

  # Call amr_search task to perform the analysis
  call run_amr_search.amr_search {
    input:
      input_fasta = input_fasta,
      samplename = samplename,
      amr_search_database = amr_search_database
  }

  output {
    File amr_search_results = amr_search.json_output
    File amr_results_csv = amr_search.output_csv
    File amr_results_pdf = amr_search.output_pdf
    String amr_search_docker = amr_search.amr_search_docker
    String amr_search_version = read_string(amr_search.output_version)
  }
}
