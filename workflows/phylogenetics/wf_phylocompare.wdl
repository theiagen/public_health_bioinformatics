version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/phylogenetic_inference/utilities/task_phylovalidate.wdl" as task_phylovalidate
import "../../tasks/phylogenetic_inference/utilities/task_cophylogeny.wdl" as task_cophylogeny
import "../../tasks/phylogenetic_inference/utilities/task_root_phylo.wdl" as task_root_phylo

workflow phylocompare {
  input {
    File tree1
    File tree2

    Boolean? validate = false

    String? outgroup
    Boolean? midpoint = false
  }
  call versioning.version_capture {
    input:
  }
  if (defined(outgroup) || midpoint == true) {
    call task_root_phylo.root_phylo as root_tree1_task {
      input:
        tree = tree1,
        outgroup = outgroup,
        midpoint = midpoint
    }
    call task_root_phylo.root_phylo as root_tree2_task {
      input:
        tree = tree2,
        outgroup = outgroup,
        midpoint = midpoint
    }
  }
  call task_cophylogeny.cophylogeny as cophylo_task {
    input:
        tree1 = select_first([root_tree1_task.rooted_tree, tree1]),
        tree2 = select_first([root_tree2_task.rooted_tree, tree2])
  }
  if (validate == true) {
    call task_phylovalidate.phylovalidate as phylovalidate_task {
      input:
          tree1 = select_first([root_tree1_task.rooted_tree, tree1]),
          tree2 = select_first([root_tree2_task.rooted_tree, tree2]),
    }
  }
  output {
    String phylocompare_phb_version = version_capture.phb_version
    String? phylovalidate_version = phylovalidate_task.phylovalidate_version
    File? phylovalidate_report = phylovalidate_task.summary_report
    String? phylovalidate_distance = phylovalidate_task.phylovalidate_distance
    String? phylovalidate_flag = phylovalidate_task.phylovalidate_flag
    File? phylovalidate_tree1_clean = phylovalidate_task.tree1_clean
    File? phylovalidate_tree2_clean = phylovalidate_task.tree2_clean
    File cophylogeny_plot_with_branch_lengths = cophylo_task.cophylogeny_with_branch_lengths
    File cophylogeny_plot = cophylo_task.cophylogeny_branch_order
    String cophylogeny_version = cophylo_task.theiaphylo_version
  }
}
