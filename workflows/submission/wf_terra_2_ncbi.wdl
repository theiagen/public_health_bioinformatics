version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/submission/task_broad_ncbi_tools.wdl" as ncbi_tools
import "../../tasks/utilities/submission/task_submission.wdl" as submission

workflow Terra_2_NCBI {
  input {
    String project_name
    String workspace_name
    String table_name
    Array[String] sample_names
    Boolean skip_biosample = false 
    File ncbi_config_js
    File? input_table # for command line testing only
    String biosample_package # used to be biosample_type
    String sra_transfer_gcp_bucket # used to be gcp_bucket_uri
    Boolean submit_to_production = false # used to be path_on_ftp_server
    String bioproject
  }
  call versioning.version_capture{
    input:
  }
  call submission.prune_table {
    input:
      project_name = project_name,
      workspace_name = workspace_name,
      table_name = table_name,
      sample_names = sample_names,
      input_table = input_table,
      biosample_type = biosample_package,
      bioproject = bioproject,
      gcp_bucket_uri = sra_transfer_gcp_bucket,
      skip_biosample = skip_biosample
  }
  if (skip_biosample == false){
    call ncbi_tools.biosample_submit_tsv_ftp_upload {
      input:
        meta_submit_tsv = prune_table.biosample_table, 
        config_js = ncbi_config_js, 
        submit_to_production = submit_to_production
    }
    call submission.add_biosample_accessions {
      input:
        generated_accessions = biosample_submit_tsv_ftp_upload.generated_accessions,
        sra_metadata = prune_table.sra_table_for_biosample,
        project_name = project_name,
        workspace_name = workspace_name,
        table_name = table_name
    }
  }
  if (select_first([add_biosample_accessions.proceed, true]) == true) {
    call ncbi_tools.sra_tsv_to_xml {
      input:
        meta_submit_tsv = select_first([add_biosample_accessions.sra_table, prune_table.sra_table]),
        config_js = ncbi_config_js,
        bioproject = bioproject,
        data_bucket_uri = sra_transfer_gcp_bucket
    }
    call ncbi_tools.ncbi_sftp_upload {
      input: 
        submission_xml = sra_tsv_to_xml.submission_xml,
        config_js = ncbi_config_js,
        submit_to_production = submit_to_production
    }
  }
  output {
    # Workflow produced files
    File sra_metadata = select_first([add_biosample_accessions.sra_table, prune_table.sra_table])
    File biosample_metadata = prune_table.biosample_table
    File excluded_samples = prune_table.excluded_samples
    File? generated_accessions = biosample_submit_tsv_ftp_upload.generated_accessions
    File? biosample_failures = biosample_submit_tsv_ftp_upload.biosample_failures
    String? biosample_status = add_biosample_accessions.biosample_status
    # NCBI produced files
    File? biosample_submission_xml = biosample_submit_tsv_ftp_upload.submission_xml
    Array[File]? biosample_report_xmls = biosample_submit_tsv_ftp_upload.report_xmls
    File? sra_submission_xml = sra_tsv_to_xml.submission_xml
    Array[File]? sra_report_xmls = ncbi_sftp_upload.reports_xmls
    # Versioning
    String Terra_2_NCBI_version = version_capture.phb_version
    String Terra_2_NCBI_analysis_date = version_capture.date
  }
}