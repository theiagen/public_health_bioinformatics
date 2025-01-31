version 1.0

import "../../tasks/gene_typing/drug_resistance/task_amr_search.wdl" as run_amr_search
import "../../tasks/utilities/data_handling/parse_amr_json.wdl" as parse_amr_json

workflow amr_search_workflow {
  input {
    File input_fasta
    String? database
    String samplename
  }

  # Call amr_search task to perform the analysis
  call run_amr_search.amr_search {
    input:
      input_fasta = input_fasta,
      samplename = samplename,
      database = database
  }

  # Call parse_amr_json task to process the output JSON
  call parse_amr_json.parse_amr_json as parse_json {
    input:
      input_json = amr_search.json_output,
      samplename = samplename,
  }

  output {
    File amr_search_results = amr_search.json_output
    File amr_results_csv = parse_json.output_csv
    File amr_results_png = parse_json.output_png
    String amr_search_docker = amr_search.amr_search_docker
    String amr_search_version = read_string(parse_json.output_version)
  }
}
