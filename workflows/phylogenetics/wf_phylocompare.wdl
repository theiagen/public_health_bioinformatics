version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/data_handling/task_phylocompare.wdl" as phylovalidate

workflow phylocompare {
  input {
    String tree1_path
    String tree2_path

    String? root_tips
    Boolean? midpoint = false
    Boolean? unrooted = false 
    Float? max_distance = 0.0
  }
  call versioning.version_capture {
    input:
  }
  call phylovalidate.phylovalidate {
    input:
        tree1_path = tree1_path,
        tree2_path = tree2_path,
        root_tips = root_tips,
        max_distance = max_distance,
        unrooted = unrooted
  }
  output {
    String phb_version = version_capture.phb_version
    String phylocompare_version = phylovalidate.phylocompare_version
    File phylocompare_report = phylovalidate.summary_report
    Float phylo_distance = phylovalidate.phylo_distance
    String validation = phylovalidate.validation
  }
}
