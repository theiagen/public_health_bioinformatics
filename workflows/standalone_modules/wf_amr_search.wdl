version 1.0

import "../../tasks/gene_typing/drug_resistance/task_amr_search.wdl" as run_amr_search
import "../../tasks/task_versioning.wdl" as versioning

workflow amr_search_workflow {
  meta {
    description: "A standalone version of Pathogenwatch's AMR profiling functionality utilizing AMRsearch tool from Pathogenwatch."
  }
  input {
    File input_fasta
    String amr_search_database
    String samplename
  }
  call run_amr_search.amr_search {
    input:
      input_fasta = input_fasta,
      samplename = samplename,
      amr_search_database = amr_search_database
  }
  call versioning.version_capture {
    input:
  }
  output {
    File amr_search_results = amr_search.amr_search_json_output
    File amr_results_csv = amr_search.amr_search_output_csv
    File amr_results_pdf = amr_search.amr_search_output_pdf
    String amr_search_all_resistances = amr_search.amr_search_all_resistances
    String amr_search_associated_resistances = amr_search.amr_search_associated_resistances
    String amr_search_docker = amr_search.amr_search_docker_image
    String amr_search_version = amr_search.amr_search_version

    # PHB Versioning
    String amr_search_wf_analysis_date = version_capture.date
    String amr_search_wf_version = version_capture.phb_version
  }
}
