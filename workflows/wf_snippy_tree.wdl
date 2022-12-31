version 1.0

import "../tasks/phylogenetic_inference/task_snippy_core.wdl" as snippy_core
import "../tasks/phylogenetic_inference/task_iqtree.wdl" as iqtree
import "../tasks/phylogenetic_inference/task_snp_dists.wdl" as snp_dists
import "../tasks/phylogenetic_inference/task_gubbins.wdl" as gubbins
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
  }
  	call versioning.version_capture{
    input:
  }
  output {
    String snippy_tree_version = version_capture.phbg_version
    String snippy_tree_snippy_version = snippy_core.snippy_version
    File snippy_tree_core_alignment = snippy_core.snippy_core_alignment
    File snippy_tree_full_alignment = snippy_core.snippy_full_alignment
    File snippy_tree_clean_full_alignment = snippy_core.snippy_full_alignment_clean
    File snippy_tree_ref = snippy_core.snippy_ref
    File snippy_tree_all_snps = snippy_core.snippy_core_tab
    File snippy_tree_snps_summary = snippy_core.snippy_txt
    File snippy_tree_vcf = snippy_core.snippy_vcf
    File? snippy_tree_iqtree = iqtree.ml_tree
    String? snippy_tree_iqtree_version = iqtree.version
    String? snippy_tree_snpdists_gubbins_version = snp_dists_gubbins.version
    File? snippy_tree_snpdists_gubbins_matrix = snp_dists_gubbins.snp_matrix
    File? snippy_tree_snpdists_gubbins_list = snp_dists_gubbins.snp_dists_molten_ordered
    String? snippy_tree_snpdists_iqtree_version = snp_dists_iqtree.version
    File? snippy_tree_snpdists_iqtree_matrix = snp_dists_iqtree.snp_matrix
    File? snippy_tree_snpdists_iqtree_list = snp_dists_iqtree.snp_dists_molten_ordered
    File? snippy_tree_gubbins_tree = gubbins.gubbins_final_tree
    File? snippy_tree_gubbins_labelled_tree = gubbins.gubbins_final_labelled_tree
    File? snippy_tree_gubbins_polymorphic_fasta = gubbins.gubbins_polymorphic_fasta
    File? snippy_tree_gubbins_recombination_gff = gubbins.gubbins_recombination_gff
    File? snippy_tree_gubbins_branch_stats = gubbins.gubbins_branch_stats
    String? snippy_tree_gubbins_version = gubbins.version
    File? snippy_tree_gubbins_timetree = gubbins.gubbins_timetree
    File? snippy_tree_gubbins_timetree_stats = gubbins.gubbins_timetree_stats
  }
}