version 1.0

import "../tasks/phylogenetic_inference/task_snippy_core.wdl" as snippy_core
import "../tasks/phylogenetic_inference/task_iqtree.wdl" as iqtree
import "../tasks/phylogenetic_inference/task_snp_dists.wdl" as snp_dists
import "../tasks/phylogenetic_inference/task_gubbins.wdl" as gubbins
import "../tasks/utilities/task_summarize_data.wdl" as data_summary
import "../tasks/task_versioning.wdl" as versioning

workflow snippy_tree_wf {
  meta {
    description: "Perform phylogenetic tree inference using iqtree (default) or snp-dist"
  }
  input {
    String tree_name
    Array[File] snippy_variants_outdir_tarball
    Array[String] samplenames
    File reference
    Boolean use_gubbins = false
    String? data_summary_terra_project
    String? data_summary_terra_workspace
    String? data_summary_terra_table
    String? data_summary_column_names # space delimited
  }
  call snippy_core.snippy_core {
    input:
      snippy_variants_outdir_tarball = snippy_variants_outdir_tarball,
      samplenames = samplenames,
      reference = reference,
      tree_name = tree_name
  }
  if (use_gubbins) {
    call gubbins.gubbins {
      input:
        alignment = snippy_core.snippy_full_alignment_clean,
        cluster_name = tree_name
    }
    call snp_dists.snp_dists as snp_dists_gubbins {
      input:
        alignment = gubbins.gubbins_polymorphic_fasta,
        cluster_name = tree_name
    }
    call snp_dists.reorder_matrix as reorder_matrix_gubbins {
      input:
        tree = gubbins.gubbins_final_labelled_tree,
        matrix = snp_dists_gubbins.snp_matrix,
        cluster_name = tree_name 
    }
  }
  if (!use_gubbins) {
    call iqtree.iqtree {
      input:
        alignment = snippy_core.snippy_full_alignment_clean,
        cluster_name = tree_name
    }
    call snp_dists.snp_dists as snp_dists_iqtree {
      input:
        alignment = snippy_core.snippy_full_alignment_clean,
        cluster_name = tree_name
    }
    call snp_dists.reorder_matrix as reorder_matrix_iqtree {
      input:
        tree = iqtree.ml_tree,
        matrix = snp_dists_iqtree.snp_matrix,
        cluster_name = tree_name 
    }
  }
  if (defined(data_summary_column_names)) {
    call data_summary.summarize_data {
      input:
        sample_names = samplenames,
        terra_project = data_summary_terra_project,
        terra_workspace = data_summary_terra_workspace,
        terra_table = data_summary_terra_table,
        column_names = data_summary_column_names
    }
  }
  call versioning.version_capture{
    input:
  }
  output {
    # version capture
    String snippy_tree_version = version_capture.phbg_version
    String snippy_tree_analysis_date = version_capture.date
    # snippy core outputs
    String snippy_tree_snippy_version = snippy_core.snippy_version
    File snippy_tree_core_alignment = snippy_core.snippy_core_alignment
    File snippy_tree_full_alignment = snippy_core.snippy_full_alignment
    File snippy_tree_clean_full_alignment = snippy_core.snippy_full_alignment_clean
    File snippy_tree_ref = snippy_core.snippy_ref
    File snippy_tree_all_snps = snippy_core.snippy_core_tab
    File snippy_tree_snps_summary = snippy_core.snippy_txt
    File snippy_tree_vcf = snippy_core.snippy_vcf
    # iqtree outputs
    String? snippy_tree_iqtree_version = iqtree.version
    # gubbins outputs
    String? snippy_tree_gubbins_version = gubbins.version
    File? snippy_tree_gubbins_labelled_tree = gubbins.gubbins_final_labelled_tree
    File? snippy_tree_gubbins_polymorphic_fasta = gubbins.gubbins_polymorphic_fasta
    File? snippy_tree_gubbins_recombination_gff = gubbins.gubbins_recombination_gff
    File? snippy_tree_gubbins_branch_stats = gubbins.gubbins_branch_stats
    File? snippy_tree_gubbins_timetree = gubbins.gubbins_timetree
    File? snippy_tree_gubbins_timetree_stats = gubbins.gubbins_timetree_stats
    # snpdists outputs
    String snippy_tree_snpdists_version = select_first([snp_dists_gubbins.version, snp_dists_iqtree.version])
    # reorder matrix outputs
    File? snippy_tree_gubbins_matrix = reorder_matrix_gubbins.ordered_midpoint_matrix
    File? snippy_tree_gubbins_tree = reorder_matrix_gubbins.midpoint_rooted_tree
    File? snippy_tree_iqtree_matrix = reorder_matrix_iqtree.ordered_midpoint_matrix
    File? snippy_tree_iqtree_tree = reorder_matrix_iqtree.midpoint_rooted_tree
    # data summary outputs
    File? snippy_tree_summarized_data = summarize_data.summarized_data
    
    }
}