version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/data_handling/task_phylocompare.wdl" as phylovalidate

workflow phylocompare {
  input {
    String tree1_path
    String tree2_path

    String? root_tips
    Boolean? unrooted = false 
    Float? lrm_max_dist = 0.0
  }
  call versioning.version_capture {
    input:
  }
  call phylovalidate.phylocompare {
    input:
        tree1_path = tree1_path,
        tree2_path = tree2_path,
        root_tips = root_tips,
        lrm_max_dist = lrm_max_dist,
        unrooted = unrooted
  }
  output {
    String phb_version = version_capture.phb_version
    String phylocompare_version = phylocompare.phylocompare_version
    File phylocompare_report = phylocompare.summary_report
    Float lrm_distance = phylocompare.lrm_distance
    String validation = phylocompare.phylovalidate
  }
}
