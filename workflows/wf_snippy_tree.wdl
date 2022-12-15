version 1.0

import "../tasks/phylogenetic_inference/task_iqtree.wdl" as iq_tree
import "../tasks/phylogenetic_inference/task_snp_dists.wdl" as snpdists

workflow snippy_tree_wf {
  meta {
    description: "Perform phylogenetic tree inference using iqtree (default) or snp-dist"
  }
  input {
    File alignment
    String cluster_name

  }
  call iq_tree.iqtree {
    input:
      alignment = alignment,
      cluster_name = cluster_name
  }
  call snpdists.snp_dists{
    input:
      alignment = alignment,
      cluster_name = cluster_name
  }
  output {
    File snippy_iqtree = iqtree.ml_tree
    String snippy_iqtree_version = iqtree.version
    String snippy_snpdists_version = snp_dists.version
    File snippy_snpdists_matrix = snp_dists.snp_matrix
    File snippy_snpdists_molten_ordered = snp_dists.snp_dists_molten_ordered
  }
}