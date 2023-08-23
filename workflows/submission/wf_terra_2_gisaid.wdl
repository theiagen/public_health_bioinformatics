version 1.0

import "../../tasks/utilities/submission/task_gisaid_cli.wdl" as gisaid_cli
import "../../tasks/task_versioning.wdl" as versioning

workflow Terra_2_GISAID {
  input {
    File concatenated_fastas
    File concatenated_metadata
    String client_id
  }
  call gisaid_cli.gisaid_upload {
    input:
      concatenated_fastas = concatenated_fastas,
      concatenated_metadata = concatenated_metadata,
      client_id = client_id
  }
  call versioning.version_capture {
  }
  output {
    String gisaid_cli_version = gisaid_upload.gisaid_cli_version
    File gisaid_logs = gisaid_upload.gisaid_logs
    File failed_uploads = gisaid_upload.failed_uploads
    String terra_2_gisaid_version = version_capture.phb_version
    String terra_2_gisaid_analysis_date = version_capture.date
  }
}