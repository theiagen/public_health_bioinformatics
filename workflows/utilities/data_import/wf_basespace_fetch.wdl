version 1.0

import "../../../tasks/task_versioning.wdl" as versioning_task
import "../../../tasks/utilities/data_import/task_basespace_api.wdl" as basespace

workflow basespace_fetch {
  input {
    String sample_name
    String basespace_sample_name
    String basespace_collection_id
    String basespace_access_token

    String? basespace_api_url
    Boolean? validate_paired_end
    Boolean? group_by_lane
  }
  call basespace.fetch_bs {
    input:
      sample_name = sample_name,
      basespace_sample_name = basespace_sample_name,
      basespace_collection_id = basespace_collection_id,
      basespace_access_token = basespace_access_token,
      basespace_api_url = basespace_api_url,
      validate_paired_end = validate_paired_end,
      group_by_lane = group_by_lane,
  }
  call versioning_task.version_capture {
    input:
  }
  output {
    String basespace_fetch_version = version_capture.phb_version
    String basespace_fetch_analysis_date = version_capture.date
    File read1 = fetch_bs.read1
    File read2 = fetch_bs.read2
    File basespace_log = fetch_bs.basespace_log
  }
}
