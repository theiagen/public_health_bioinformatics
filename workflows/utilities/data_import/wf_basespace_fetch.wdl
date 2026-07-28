version 1.0

import "../../../tasks/task_versioning.wdl" as versioning_task
import "../../../tasks/utilities/data_import/task_basespace_api.wdl" as basespace

workflow basespace_fetch {
  input {
    String sample_name
    String basespace_sample_name
    String basespace_collection_id
    String access_token

    Boolean? validate_paired_end
  }
  call basespace.fetch_bs {
    input:
      sample_name = sample_name,
      basespace_sample_name = basespace_sample_name,
      basespace_collection_id = basespace_collection_id,
      access_token = access_token,
      validate_paired_end = validate_paired_end,
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
