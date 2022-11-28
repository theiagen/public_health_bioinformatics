version 1.0

import "../tasks/phylogenetic_inference/task_ksnp3.wdl" as ksnp3
import "../tasks/phylogenetic_inference/task_snp_dists.wdl" as snp_dists
import "../tasks/task_versioning.wdl" as versioning

workflow ksnp3_workflow {
  input {
    Array[File] assembly_fasta
    Array[String] samplename
    String cluster_name
	}
	call ksnp3.ksnp3 as ksnp3_task {
		input:
			assembly_fasta = assembly_fasta,
      samplename = samplename,
      cluster_name = cluster_name
  }
  call snp_dists.snp_dists as core_snp_dists {
    input:
      cluster_name = cluster_name,
      alignment = ksnp3_task.ksnp3_core_matrix
  }
  call snp_dists.snp_dists as pan_snp_dists {
    input:
      cluster_name = cluster_name,
      alignment = ksnp3_task.ksnp3_pan_matrix
  }
  call versioning.version_capture{
    input:
  }
  output {
    # Version Capture
    String ksnp3_wf_version = version_capture.phbg_version
    String ksnp3_wf_analysis_date = version_capture.date
    # ksnp3_outputs
    String ksnp3_snp_dists_version = pan_snp_dists.version
    File ksnp3_core_snp_matrix = core_snp_dists.snp_matrix
    File ksnp3_core_tree = ksnp3_task.ksnp3_core_tree
    File ksnp3_core_vcf = ksnp3_task.ksnp3_core_vcf
    File ksnp3_pan_snp_matrix = pan_snp_dists.snp_matrix
    File ksnp3_pan_parsimony_tree = ksnp3_task.ksnp3_pan_parsimony_tree
    File? ksnp3_ml_tree = ksnp3_task.ksnp3_ml_tree
    File? ksnp3_nj_tree = ksnp3_task.ksnp3_nj_tree
    String ksnp3_docker = ksnp3_task.ksnp3_docker_image
  }
}
