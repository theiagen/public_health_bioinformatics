version 1.0

import "../../../tasks/utilities/task_basespace_cli.wdl" as basespace
import "../../../tasks/task_versioning.wdl" as versioning_task

workflow basespace_fetch {
  input {
    String sample_name
    String basespace_sample_name
    String? basespace_sample_id
    String basespace_collection_id
    String api_server
    String access_token
  }
  call basespace.fetch_bs {
    input:
      sample_name = sample_name,
      basespace_sample_id = basespace_sample_id,
      basespace_sample_name = basespace_sample_name,
      basespace_collection_id = basespace_collection_id,
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