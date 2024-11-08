version 1.0

import "../../tasks/phylogenetic_inference/task_mashtree.wdl" as mashtree
import "../../tasks/phylogenetic_inference/utilities/task_reorder_matrix.wdl" as reorder_matrix_task
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/data_handling/task_summarize_data.wdl" as data_summary

workflow mashtree_fasta {
  input {
    Array[File] assembly_fasta
    String cluster_name
    Array[String]? sample_names
    String? data_summary_terra_project
    String? data_summary_terra_workspace
    String? data_summary_terra_table
    String? data_summary_column_names
    Boolean midpoint_root_tree = true
    Boolean phandango_coloring = false
  }
  String cluster_name_updated = sub(cluster_name, " ", "_")
  call mashtree.mashtree_fasta as mashtree_task {
    input:
      assembly_fasta = assembly_fasta,
      cluster_name = cluster_name_updated
    }
  call reorder_matrix_task.reorder_matrix {
    input:
      input_tree = mashtree_task.mashtree_tree,
      matrix = mashtree_task.mashtree_matrix,
      cluster_name = cluster_name_updated,
      midpoint_root_tree = midpoint_root_tree,
      phandango_coloring = phandango_coloring
  }
  if (defined(data_summary_column_names)) {
    call data_summary.summarize_data {
      input:
        sample_names = sample_names,
        terra_project = data_summary_terra_project,
        terra_workspace = data_summary_terra_workspace,
        terra_table = data_summary_terra_table,
        column_names = data_summary_column_names,
        output_prefix = cluster_name_updated,
        phandango_coloring = phandango_coloring
    }
  } 
  call versioning.version_capture {
    input:
  }
  output {
    # Versioning
    String mashtree_wf_version = version_capture.phb_version
    String mashtree_wf_analysis_date = version_capture.date
    # Masthree Out
    File mashtree_matrix = reorder_matrix.ordered_matrix
    File mashtree_tree = reorder_matrix.tree
    String mashtree_version = mashtree_task.version
    String mashtree_docker = mashtree_task.mashtree_docker
    # Data Summary Out
    File? mashtree_summarized_data = summarize_data.summarized_data
    File? mashtree_filtered_metadata = summarize_data.filtered_metadata
  }
}