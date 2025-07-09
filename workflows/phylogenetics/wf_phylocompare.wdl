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
    call task_root_phylo.root_phylo as root_tree1_task {
      input:
        tree = tree1,
        outgroups = outgroups,
        midpoint = midpoint
    }
    call task_root_phylo.root_phylo as root_tree2_task {
      input:
        tree = tree2,
        outgroups = outgroups,
        midpoint = midpoint
    }
  }
  call task_phylocompare.phylovalidate as phylovalidate_task {
    input:
        tree1 = select_first([root_tree1_task.rooted_tree, tree1]),
        tree2 = select_first([root_tree2_task.rooted_tree, tree2]),
        max_distance = max_distance
  }
  output {
    String phylocompare_phb_version = version_capture.phb_version
    String phylocompare_version = phylovalidate.phylocompare_version
    File phylocompare_report = phylovalidate.summary_report
    String phylocompare_distance = phylovalidate.phylovalidate_distance
    String phylocompare_validation = phylovalidate.phylovalidate_validation
    String phylocompare_flag = phylovalidate.phylovalidate_flag
    File phylocompare_tree1_clean = phylovalidate.tree1_clean
    File phylocompare_tree2_clean = phylovalidate.tree2_clean
  }
}
