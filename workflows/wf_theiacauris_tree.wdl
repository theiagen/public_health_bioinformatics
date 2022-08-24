version 1.0

import "../tasks/phylogenetic_inference/task_theiacauris_tree.wdl" as mashtree
import "../tasks/task_versioning.wdl" as versioning

workflow theiacauris_tree {
    input {
    File assembly_fasta
    String cluster_name = "CladeTyper_Tree"
    }
    call mashtree.theiacauris_mashtree_fasta as mashtree_task {
        input:
        assembly_fasta = assembly_fasta,
        cluster_name = cluster_name
    }
    call versioning.version_capture {
    input:
  }
  output {
        # Versioning
    String theiacauris_tree_wf_version = version_capture.phbg_version
    String theiacauris_tree_wf_analysis_date = version_capture.date
        # Masthree Out
    File theiacauris_tree_matrix = mashtree_task.mashtree_matrix
    File theiacauris_tree_tree = mashtree_task.mashtree_tree
    String theiacauris_tree_version = mashtree_task.version
    }
}