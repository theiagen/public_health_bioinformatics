version 1.0

import "../../tasks/variant_calling/task_clair3.wdl" as clair3
import "../../tasks/alignment/task_minimap2.wdl" as minimap2
import "../../tasks/task_versioning.wdl" as versioning

workflow clair3phy_workflow {
  input {
    Array[File] read1
    Array[File] read2
    File reference
    String cluster_name
    Boolean midpoint_root_tree = true
	}
  call versioning.version_capture {
    input:
  }
  output {
    # Version Capture
    String clair3phy_wf_version = version_capture.phb_version
    String clair3phy_wf_analysis_date = version_capture.date
  }
}