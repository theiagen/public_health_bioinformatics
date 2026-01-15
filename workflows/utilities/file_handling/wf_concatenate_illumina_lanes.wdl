version 1.0

import "../../../tasks/utilities/file_handling/task_cat_lanes.wdl" as concatenate_lanes
import "../../../tasks/task_versioning.wdl" as versioning

workflow concatenate_illumina_lanes {
  input {
    String samplename
    
    File read1_lane1
    File read1_lane2
    File? read1_lane3
    File? read1_lane4
    
    File? read2_lane1
    File? read2_lane2
    File? read2_lane3
    File? read2_lane4
  }
  call concatenate_lanes.cat_lanes {
    input:
      samplename = samplename,
      read1_lane1 = read1_lane1,
      read2_lane1 = read2_lane1,
      read1_lane2 = read1_lane2,
      read2_lane2 = read2_lane2,
      read1_lane3 = read1_lane3,
      read2_lane3 = read2_lane3,
      read1_lane4 = read1_lane4,
      read2_lane4 = read2_lane4
  }
  call versioning.version_capture {
    input:
  }
  output {
    String concatenate_illumina_lanes_version = version_capture.phb_version
    String concatenate_illumina_lanes_analysis_date = version_capture.date

    File read1_concatenated = cat_lanes.read1_concatenated
    File? read2_concatenated = cat_lanes.read2_concatenated
  }
}