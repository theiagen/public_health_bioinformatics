version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/phylogenetic_inference/utilities/task_phylocompare.wdl" as task_phylocompare
import "../../tasks/phylogenetic_inference/utilities/task_root_phylo.wdl" as task_root_phylo

workflow phylocompare {
  input {
    String tree1
    String tree2

    String? outgroups
    Boolean? midpoint = false
    Float? max_distance = 0.0
  }
  call versioning.version_capture {
    input:
  }
  if (defined(outgroups) || midpoint == true) {
    call task_root_phylo.root_phylo as root1 {
      input:
        tree = tree1,
        outgroups = outgroups,
        midpoint = midpoint
    }
    call task_root_phylo.root_phylo as root2 {
      input:
        tree = tree2,
        outgroups = outgroups,
        midpoint = midpoint
    }
  }
  call task_phylocompare.phylovalidate {
    input:
        tree1 = select_first([root1.rooted_tree, tree1]),
        tree2 = select_first([root2.rooted_tree, tree2]),
        max_distance = max_distance
  }
  output {
    String phb_version = version_capture.phb_version
    String phylocompare_version = phylovalidate.phylocompare_version
    File phylocompare_report = phylovalidate.summary_report
    String phylo_distance = phylovalidate.phylo_distance
    String phylo_validation = phylovalidate.phylo_validation
    String phylo_flag = phylovalidate.phylo_flag
    File tree1_final = phylovalidate.tree1_clean
    File tree2_final = phylovalidate.tree2_clean
  }
}
