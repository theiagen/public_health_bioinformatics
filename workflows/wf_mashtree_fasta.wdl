version 1.0

import "../tasks/phylogenetic_inference/task_mashtree.wdl" as mashtree
import "../tasks/task_versioning.wdl" as versioning

workflow mashtree_fasta {
	input {
		Array[File] assembly_fasta
    String cluster_name
	}
	call mashtree.mashtree_fasta as mashtree_task {
		input:
			assembly_fasta = assembly_fasta,
      cluster_name = cluster_name
	}
	call versioning.version_capture{
    input:
  }
  output {
		# Versioning
    String mashtree_wf_version = version_capture.phbg_version
    String mashtree_wf_analysis_date = version_capture.date
		# Masthree Out
    File mashtree_matrix = mashtree_task.mashtree_matrix
		File mashtree_tree = mashtree_task.mashtree_tree
		String mashtree_version = mashtree_task.version
    }
}