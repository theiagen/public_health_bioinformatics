version 1.0

import "../tasks/phylogenetic_inference/task_pirate.wdl" as pirate
import "../tasks/phylogenetic_inference/task_iqtree.wdl" as iqtree
import "../tasks/phylogenetic_inference/task_snp_dists.wdl" as snp_dists
import "../tasks/task_versioning.wdl" as versioning

workflow pangenome_snp_workflow {
  input {
    Array[File] prokka_gff
    String cluster_name
    Boolean align = true
    Boolean core = true
  }
  call pirate.pirate as pirate {
    input:
      prokka_gff = prokka_gff,
      cluster_name = cluster_name,
      align = align
  }
  if (align == true) {
    call iqtree.iqtree as pan_iqtree {
      input:
        alignment = pirate.pirate_pangenome_alignment_fasta,
        cluster_name = cluster_name
    }
    call snp_dists.snp_dists as pan_snp_dists {
      input:
        alignment = pirate.pirate_pangenome_alignment_fasta,
        cluster_name = cluster_name
    }
    if (core == true) {
      call iqtree.iqtree as core_iqtree {
        input:
          alignment = pirate.pirate_core_alignment_fasta,
          cluster_name = cluster_name
      }
      call snp_dists.snp_dists as core_snp_dists {
        input:
          alignment = pirate.pirate_core_alignment_fasta,
          cluster_name = cluster_name
      }
    }
  }
  call versioning.version_capture{
    input:
  }
  output {
    # Version Capture
    String pirate_wf_version = version_capture.phbg_version
    String pirate_wf_analysis_date = version_capture.date
    # pirate_outputs
    File pirate_pangenome_summary = pirate.pirate_pangenome_summary
    File pirate_gene_families_ordered = pirate.pirate_gene_families_ordered
    #File pirate_unique_alleles = pirate.pirate_unique_alleles
    #File pirate_binary_fasta = pirate.pirate_binary_fasta
    #File pirate_binary_tree = pirate.pirate_binary_tree
    #File pirate_pangenome_gfa = pirate.pirate_pangenome_gfa
    #File? pirate_pangenome_alignment_fasta = pirate.pirate_pangenome_alignment_fasta
    #File? pirate_pangenome_alignment_gff = pirate.pirate_pangenome_alignment_gff
    #File? pirate_core_alignment_fasta = pirate.pirate_core_alignment_fasta
    #File? pirate_core_alignment_gff = pirate.pirate_core_alignment_gff
    String pirate_docker_image = pirate.pirate_docker_image
    # snp_dists outputs
    String? pirate_snps_dists_version = pan_snp_dists.version
    File? pirate_core_snp_matrix = core_snp_dists.snp_matrix
    File? pirate_pan_snp_matrix = pan_snp_dists.snp_matrix
    # iqtree outputs
    String? pirate_iqtree_version = pan_iqtree.version
    File? pirate_iqtree_core_tree = core_iqtree.ml_tree
    File? pirate_iqtree_pan_tree = pan_iqtree.ml_tree
  }
}