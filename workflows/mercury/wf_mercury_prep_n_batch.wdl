version 1.0

import "../tasks/task_mercury_file_wrangling.wdl" as submission
import "../tasks/task_versioning.wdl" as versioning

workflow mercury_prep_n_batch {
  input {
    String table_name
    String workspace_name
    String project_name
    Array[String] sample_names
    String organism = "sars-cov-2"
    String output_name
    String gcp_bucket_uri
    File? input_table
    Int vadr_alert_limit = 0 # only for SC2
    File? authors_sbt # only for mpox
    Boolean skip_county = false
    Boolean skip_ncbi = false
  }
  call submission.sm_metadata_wrangling {
    input:
      table_name = table_name,
      workspace_name = workspace_name,
      project_name = project_name,
      sample_names = sample_names,
      organism = organism, 
      output_name = output_name,
      gcp_bucket_uri = gcp_bucket_uri,
      input_table = input_table,
      vadr_alert_limit = vadr_alert_limit,
      skip_county = skip_county,
      skip_ncbi = skip_ncbi
  }
  if (organism == "sars-cov-2" && skip_ncbi == false) {
    call submission.trim_genbank_fastas {
      input:
        genbank_untrimmed_fasta = select_first([sm_metadata_wrangling.genbank_untrimmed_fasta]),
        output_name = output_name
    }
  }
  if (organism == "mpox" && skip_ncbi == false) {
    call submission.table2asn {
      input:
        authors_sbt = select_first([authors_sbt]),
        bankit_fasta = select_first([sm_metadata_wrangling.bankit_fasta]),
        bankit_metadata = select_first([sm_metadata_wrangling.bankit_metadata]),
        output_name = output_name
    }
  }
  call versioning.version_capture {
  }
  output {
    File excluded_samples = sm_metadata_wrangling.excluded_samples
    File? biosample_metadata = sm_metadata_wrangling.biosample_metadata
    File? sra_metadata = sm_metadata_wrangling.sra_metadata
    File? genbank_metadata = sm_metadata_wrangling.genbank_metadata
    File? genbank_fasta = trim_genbank_fastas.genbank_fasta
    File? bankit_sqn_to_email = table2asn.sqn_file
    File gisaid_metadata = sm_metadata_wrangling.gisaid_metadata
    File gisaid_fasta = sm_metadata_wrangling.gisaid_fasta
    String mercury_prep_n_batch_analysis_date = version_capture.date
    String phvg_version = version_capture.phvg_version
  }
}