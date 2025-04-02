version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/data_handling/task_phylocompare.wdl" as task_phylocompare

workflow phylocompare {
  input {
    String tree1_path
    String tree2_path

    String? outgroups
    Boolean? midpoint = false
    Float? max_distance = 0.0
  }
  call versioning.version_capture {
    input:
  }
  call task_phylocompare.phylovalidate {
    input:
        tree1_path = tree1_path,
        tree2_path = tree2_path,
        outgroups = outgroups,
        midpoint = midpoint,
        max_distance = max_distance
  }
  output {
    String phb_version = version_capture.phb_version
    String phylocompare_version = phylovalidate.phylocompare_version
    File phylocompare_report = phylovalidate.summary_report
    String phylo_distance = phylovalidate.phylo_distance
    String phylo_validation = phylovalidate.phylo_validation
    String phylo_flag = phylovalidate.phylo_flag
  }
}
