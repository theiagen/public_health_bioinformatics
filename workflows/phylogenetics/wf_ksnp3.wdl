version 1.0

import "../../tasks/phylogenetic_inference/task_ksnp3.wdl" as ksnp3
import "../../tasks/phylogenetic_inference/utilities/task_ksnp_shared_snps.wdl" as ksnp_shared_snps
import "../../tasks/phylogenetic_inference/utilities/task_reorder_matrix.wdl" as reorder_matrix
import "../../tasks/phylogenetic_inference/utilities/task_snp_dists.wdl" as snp_dists
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/data_handling/task_summarize_data.wdl" as data_summary

workflow ksnp3_workflow {
  input {
    Array[File] assembly_fasta
    Array[String] samplename
    String cluster_name
    String cluster_name_updated = sub(cluster_name, " ", "_")
    String? data_summary_terra_project
    String? data_summary_terra_workspace
    String? data_summary_terra_table
    String? data_summary_column_names # string of comma delimited column names
    Boolean midpoint_root_tree = true
    Boolean phandango_coloring = false
	}
  call ksnp3.ksnp3 as ksnp3_task {
    input:
      assembly_fasta = assembly_fasta,
      samplename = samplename,
      cluster_name = cluster_name_updated
  }
  if (ksnp3_task.skip_core_snp_dists == "The core SNP matrix was produced") {
    call snp_dists.snp_dists as core_snp_dists {
      input:
        cluster_name = cluster_name_updated,
        alignment = ksnp3_task.ksnp3_core_matrix
    }
    call reorder_matrix.reorder_matrix as core_reorder_matrix {
      input:
        input_tree = ksnp3_task.ksnp3_core_tree,
        matrix = core_snp_dists.snp_matrix,
        cluster_name = cluster_name_updated + "_core",
        midpoint_root_tree = midpoint_root_tree,
        phandango_coloring = phandango_coloring
    }
    call ksnp_shared_snps.ksnp_shared_snps as core_ksnp3_shared_snps_task {
      input:
        ksnp_vcf_ref_genome = ksnp3_task.ksnp3_vcf_ref_genome,
        samplename = samplename,
        cluster_name = cluster_name_updated
    }
  }
  call snp_dists.snp_dists as pan_snp_dists {
    input:
      cluster_name = cluster_name_updated,
      alignment = ksnp3_task.ksnp3_pan_matrix
  }
  call reorder_matrix.reorder_matrix as pan_reorder_matrix {
    input:
      input_tree = ksnp3_task.ksnp3_pan_parsimony_tree,
      matrix = pan_snp_dists.snp_matrix,
      cluster_name = cluster_name_updated + "_pan",
      midpoint_root_tree = midpoint_root_tree,
      phandango_coloring = phandango_coloring
  }
  if (defined(data_summary_column_names)) {
    call data_summary.summarize_data {
      input:
        sample_names = samplename,
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
    # Version Capture
    String ksnp3_wf_version = version_capture.phb_version
    String ksnp3_wf_analysis_date = version_capture.date
    String ksnp3_docker = ksnp3_task.ksnp3_docker_image
    # ksnp3_outputs
    String ksnp3_snp_dists_version = pan_snp_dists.snp_dists_version
    File ksnp3_vcf_ref_genome = ksnp3_task.ksnp3_vcf_ref_genome
    File ksnp3_vcf_snps_not_in_ref = ksnp3_task.ksnp3_vcf_snps_not_in_ref
    String ksnp3_vcf_ref_samplename = ksnp3_task.ksnp3_vcf_ref_samplename
    String ksnp3_core_snp_matrix_status = ksnp3_task.skip_core_snp_dists
    File ksnp3_snps = ksnp3_task.ksnp3_snps_all
    String ksnp3_number_snps = ksnp3_task.ksnp3_number_snps
    String ksnp3_number_core_snps = ksnp3_task.ksnp3_number_core_snps
    # ordered matrixes and reordered trees
    File? ksnp3_core_snp_matrix = core_reorder_matrix.ordered_matrix
    File? ksnp3_core_tree = core_reorder_matrix.tree
    File ksnp3_pan_snp_matrix = pan_reorder_matrix.ordered_matrix
    File ksnp3_pan_tree = pan_reorder_matrix.tree
    # optional tree outputs
    File? ksnp3_ml_tree = ksnp3_task.ksnp3_ml_tree
    File? ksnp3_nj_tree = ksnp3_task.ksnp3_nj_tree
    # data summary output 
    File? ksnp3_summarized_data = summarize_data.summarized_data
    File? ksnp3_filtered_metadata = summarize_data.filtered_metadata
    # ksnp3_shared_snps outputs
    File? ksnp3_core_snp_table = core_ksnp3_shared_snps_task.core_snp_table
  }
}