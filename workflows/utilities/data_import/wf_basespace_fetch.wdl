version 1.0

import "../../../tasks/task_versioning.wdl" as versioning_task
import "../../../tasks/utilities/data_import/task_basespace_cli.wdl" as basespace

workflow basespace_fetch {
  input {
    String basespace_sample_name
    String? basespace_run_id
    String? basespace_project_id
    String api_server
    String access_token
  }
  call basespace.fetch_bs {
    input:
      basespace_sample_name = basespace_sample_name,
      basespace_run_id = basespace_run_id,
      basespace_project_id = basespace_project_id,
      api_server = api_server,
      access_token = access_token
  }
  call versioning_task.version_capture {
    input:
  }
  output {
    String basespace_fetch_version = version_capture.phb_version
    String basespace_fetch_analysis_date = version_capture.date
    File read1 = fetch_bs.read1
    File? read2 = fetch_bs.read2
  }
}