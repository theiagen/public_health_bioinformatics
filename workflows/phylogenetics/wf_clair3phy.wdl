version 1.0

import "../../tasks/variant_calling/task_clair3.wdl" as clair3_task
import "../../tasks/alignment/task_minimap2.wdl" as minimap2_task
import "../../tasks/utilities/data_handling/task_parse_mapping.wdl" as parse_mapping_task
import "../../tasks/task_versioning.wdl" as versioning

workflow clair3phy_workflow {
  input {
    Array[File] read1
    Array[String] samplenames
    File reference
    String cluster_name
	}
  call versioning.version_capture {
    input:
  }
  call minimap2_task.minimap2_merge {
    input:
      query1 = read1,
      reference = reference,
      custer_name = cluster_name,
      mode = "map-ont",
      output_sam = true,
      additional_options = "-L --cs --MD"
  }
  call parse_mapping_task.sam_to_sorted_bam {
    input:
      sam = minimap2_merge.minimap2_out,
      samplename = cluster_name
  }
  call clair3_task.clair3 {
    input:
      alignment = sam_to_sorted_bam.bam,
      alignment_index = sam_to_sorted_bam.bai,
      reference = reference,
      cluster_name = cluster_name
  }
  output {
    # Version Capture
    String clair3phy_wf_version = version_capture.phb_version
    String clair3phy_wf_analysis_date = version_capture.date
  }
}