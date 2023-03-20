version 1.0

import "../../tasks/phylogenetic_inference/task_snippy_core.wdl" as snippy_core_task
import "../../tasks/phylogenetic_inference/task_snp_sites.wdl" as snp_sites_task
import "../../tasks/phylogenetic_inference/task_iqtree.wdl" as iqtree_task
import "../../tasks/phylogenetic_inference/task_snp_dists.wdl" as snp_dists
import "../../tasks/phylogenetic_inference/task_reorder_matrix.wdl" as reorder_matrix
import "../../tasks/phylogenetic_inference/task_gubbins.wdl" as gubbins_task
import "../../tasks/utilities/task_summarize_data.wdl" as data_summary
import "../../tasks/task_versioning.wdl" as versioning

workflow snippy_tree_wf {
  meta {
    description: "Perform phylogenetic tree inference using iqtree (default) or snp-dist"
  }
  input {
    String tree_name
    Array[File] snippy_variants_outdir_tarball
    Array[String] samplenames
    File reference
    Boolean use_gubbins = true
    Boolean core_genome = true
    String? data_summary_terra_project
    String? data_summary_terra_workspace
    String? data_summary_terra_table
    String? data_summary_column_names # comma delimited
  }
  call snippy_core_task.snippy_core {
    input:
      snippy_variants_outdir_tarball = snippy_variants_outdir_tarball,
      samplenames = samplenames,
      reference = reference,
      tree_name = tree_name
  }
  if (use_gubbins) {
    call gubbins_task.gubbins {
      input:
        alignment = snippy_core.snippy_full_alignment_clean,
        cluster_name = tree_name
    }
    if (core_genome) {
      call snp_sites_task.snp_sites as snp_sites_gubbins{
        input:
          msa_fasta = gubbins.gubbins_polymorphic_fasta,
          output_name = tree_name,
          output_vcf = false,
          output_multifasta = true
      }
    }
  }
  if (!use_gubbins) {
    if (core_genome) {
      call snp_sites_task.snp_sites as snp_sites_no_gubbins {
        input:
          msa_fasta = snippy_core.snippy_full_alignment_clean,
          output_name = tree_name,
          output_vcf = false,
          output_multifasta = true
      }
    }
  }
  call iqtree_task.iqtree {
    input:
      alignment = select_first([snp_sites_gubbins.snp_sites_multifasta, gubbins.gubbins_polymorphic_fasta, snp_sites_no_gubbins.snp_sites_multifasta, snippy_core.snippy_full_alignment_clean]),
      cluster_name = tree_name
  }
  call snp_dists.snp_dists as snp_dists_iqtree {
    input:
      alignment = select_first([snp_sites_gubbins.snp_sites_multifasta, gubbins.gubbins_polymorphic_fasta, snp_sites_no_gubbins.snp_sites_multifasta, snippy_core.snippy_full_alignment_clean]),
      cluster_name = tree_name
  }
  call reorder_matrix.reorder_matrix as reorder_matrix_iqtree {
    input:
      input_tree = iqtree.ml_tree,
      matrix = snp_dists_iqtree.snp_matrix,
      cluster_name = tree_name 
  }
  if (defined(data_summary_column_names)) {
    call data_summary.summarize_data {
      input:
        sample_names = samplenames,
        terra_project = data_summary_terra_project,
        terra_workspace = data_summary_terra_workspace,
        terra_table = data_summary_terra_table,
        column_names = data_summary_column_names,
        output_prefix = tree_name
    }
  }
  call versioning.version_capture{
    input:
  }
  output {
    # version capture
    String snippy_tree_version = version_capture.phb_version
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
    String snippy_tree_snpdists_version = snp_dists_iqtree.version
    # reorder matrix outputs
    File? snippy_tree_gubbins_matrix = reorder_matrix_iqtree.ordered_matrix
    File? snippy_tree_gubbins_tree = reorder_matrix_iqtree.tree
    File? snippy_tree_iqtree_matrix = reorder_matrix_iqtree.ordered_matrix
    File? snippy_tree_iqtree_tree = reorder_matrix_iqtree.tree
    # data summary outputs
    File? snippy_tree_summarized_data = summarize_data.summarized_data
    }
}