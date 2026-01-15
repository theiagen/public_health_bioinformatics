version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/data_export/task_download_terra_table.wdl" as download_table
import "../../tasks/utilities/submission/task_mercury.wdl" as mercury_task
import "../../tasks/utilities/submission/task_mercury_utilities.wdl" as submission_utilities

workflow mercury_prep_n_batch {
  input {
    String terra_table_name
    String terra_workspace_name
    String terra_project_name
    Array[String] sample_names
    String organism = "sars-cov-2"
    String output_name = "mercury"
    String gcp_bucket_uri
    File? authors_sbt # only for mpox
    Boolean skip_ncbi = false
  }
  String output_name_updated = sub(output_name, " ", "_")
  call download_table.download_terra_table {
    input:
      terra_table_name = terra_table_name,
      terra_workspace_name = terra_workspace_name,
      terra_project_name = terra_project_name
  }
  call mercury_task.mercury {
    input:
      samplenames = sample_names,
      gcp_bucket_uri = gcp_bucket_uri,
      terra_project_name = terra_project_name,
      terra_workspace_name = terra_workspace_name,
      data_table = download_terra_table.terra_table,
      table_name = terra_table_name,
      organism = organism,
      output_name = output_name_updated,
      skip_ncbi = skip_ncbi
  }
  if (mercury.organism_name == "sars-cov-2" && skip_ncbi == false) {
    call submission_utilities.trim_genbank_fastas {
      input:
        genbank_untrimmed_fasta = select_first([mercury.genbank_fasta]),
        output_name = output_name_updated
    }
  }
  if (mercury.organism_name == "mpox" && skip_ncbi == false) {
    call submission_utilities.table2asn {
      input:
        authors_sbt = select_first([authors_sbt]),
        bankit_fasta = select_first([mercury.bankit_fasta]),
        bankit_metadata = select_first([mercury.bankit_metadata]),
        output_name = output_name_updated
    }
  }
  call versioning.version_capture {
  }
  output {
    File? bankit_fasta = mercury.bankit_fasta
    File? bankit_metadata = mercury.bankit_metadata
    File? excluded_samples = mercury.excluded_samples
    File? biosample_metadata = mercury.biosample_metadata
    File? sra_metadata = mercury.sra_metadata
    File? genbank_metadata = mercury.genbank_metadata
    File? genbank_fasta = trim_genbank_fastas.genbank_fasta
    File? bankit_sqn_to_email = table2asn.sqn_file
    File? gisaid_metadata = mercury.gisaid_metadata
    File? gisaid_fasta = mercury.gisaid_fasta
    String mercury_script_version = mercury.mercury_version
    String mercury_prep_n_batch_analysis_date = version_capture.date
    String mercury_prep_n_batch_version = version_capture.phb_version
  }
}