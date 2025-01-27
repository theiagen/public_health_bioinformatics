version 1.0

import "../../tasks/gene_typing/drug_resistance/amr_search.wdl" as amr_search
import "../../tasks/utilities/data_handling/parse_amr_json.wdl" as parse_amr_json

workflow amr_search_workflow {
  input {
    File input_fasta
    String database
    String samplename
  }

  # Call amr_search task to perform the analysis
  call amr_search.amr_search {
    input:
      input_fasta = input_fasta,
      samplename = samplename,
      database = database
  }

  # Call parse_amr_json task to process the output JSON
  call parse_amr_json.parse_amr_json {
    input:
      input_json = amr_search.json_output,
      output_csv_name = samplename + "_amr_results.csv"
  }
  output {
    File amr_search_results = amr_search.json_output
    File amr_results_csv = parse_amr_json.output_csv
    File amr_results_png = parse_amr_json.output_png
  }
}
