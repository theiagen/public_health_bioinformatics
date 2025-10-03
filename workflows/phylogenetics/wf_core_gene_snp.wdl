version 1.0

import "../../tasks/phylogenetic_inference/task_iqtree.wdl" as iqtree
import "../../tasks/phylogenetic_inference/task_pirate.wdl" as pirate_task
import "../../tasks/phylogenetic_inference/utilities/task_reorder_matrix.wdl" as reorder_matrix
import "../../tasks/phylogenetic_inference/utilities/task_snp_dists.wdl" as snp_dists
import "../../tasks/phylogenetic_inference/utilities/task_snp_sites.wdl" as snp_sites
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/data_handling/task_summarize_data.wdl" as data_summary


workflow core_gene_snp_workflow {
  input {
    Array[File] gff3
    String cluster_name
    # if align = true, the pirate task will produce core and pangenome alignments for the sample set,
    # otherwise, pirate will only produce a pangenome summary
    Boolean align = true
    # use core_tree = true to produce a phylogenetic tree and snp distance matrix from the core genome alignment
    Boolean core_tree = true
    # use pan_tree = true to produce a phylogenetic tree and snp distance matrix from the pangenome alignment
    Boolean pan_tree = false
    # call snp_sites from core genome alignment
    Boolean snp_sites = true
    # data summary input variables
    Array[String]? sample_names
    String? data_summary_terra_project
    String? data_summary_terra_workspace
    String? data_summary_terra_table
    String? data_summary_column_names
    Boolean phandango_coloring = false

    Boolean midpoint_root_tree = true
  }
  String cluster_name_updated = sub(cluster_name, " ", "_")
  call pirate_task.pirate {
    input:
      gff3 = gff3,
      cluster_name = cluster_name_updated,
      align = align
  }
  if (align) {
    if (core_tree) {
      if (snp_sites) {
        call snp_sites.snp_sites {
          input:
            msa_fasta = select_first([pirate.pirate_core_alignment_fasta]),
            output_name = cluster_name_updated + "_core",
            allow_wildcard_bases = false,
            output_vcf = false,
            output_phylip = false,
            output_multifasta = true,
            output_pseudo_ref = false,
            output_monomorphic = false
        }
      }
      call iqtree.iqtree as core_iqtree {
        input:
          alignment = select_first([snp_sites.snp_sites_multifasta, pirate.pirate_core_alignment_fasta]),
          cluster_name = cluster_name_updated
      }
      call snp_dists.snp_dists as core_snp_dists {
        input:
          alignment = select_first([snp_sites.snp_sites_multifasta, pirate.pirate_core_alignment_fasta]),
          cluster_name = cluster_name_updated
      }
      call reorder_matrix.reorder_matrix as core_reorder_matrix {
        input:
          input_tree = core_iqtree.ml_tree,
          matrix = core_snp_dists.snp_matrix,
          cluster_name = cluster_name_updated + "_core",
          midpoint_root_tree = midpoint_root_tree,
          phandango_coloring = phandango_coloring
      }
    }
    if (pan_tree) {
      call iqtree.iqtree as pan_iqtree {
        input:
          alignment = select_first([pirate.pirate_pangenome_alignment_fasta]),
          cluster_name = cluster_name_updated
      }
      call snp_dists.snp_dists as pan_snp_dists {
        input:
          alignment = select_first([pirate.pirate_pangenome_alignment_fasta]),
          cluster_name = cluster_name_updated
      }
      call reorder_matrix.reorder_matrix as pan_reorder_matrix {
        input:
          input_tree = pan_iqtree.ml_tree,
          matrix = pan_snp_dists.snp_matrix,
          cluster_name = cluster_name_updated + "_pan",
          midpoint_root_tree = midpoint_root_tree,
          phandango_coloring = phandango_coloring
      }
    }
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
    # Version Capture
    String core_gene_snp_wf_version = version_capture.phb_version
    String core_gene_snp_wf_analysis_date = version_capture.date
    # pirate_outputs
    File pirate_pangenome_summary = pirate.pirate_pangenome_summary
    File pirate_gene_families_ordered = pirate.pirate_gene_families_ordered
    File? pirate_core_alignment_fasta = pirate.pirate_core_alignment_fasta
    File? pirate_core_alignment_gff = pirate.pirate_core_alignment_gff
    File? pirate_pan_alignment_fasta = pirate.pirate_pangenome_alignment_fasta
    File? pirate_pan_alignment_gff = pirate.pirate_pangenome_alignment_gff
    File? pirate_presence_absence_csv = pirate.pirate_presence_absence_csv
    String pirate_docker_image = pirate.pirate_docker_image
    # snp_sites outputs
    File? snp_sites_multifasta = snp_sites.snp_sites_multifasta
    String? snp_sites_version = snp_sites.snp_sites_version
    String? snp_sites_docker = snp_sites.snp_sites_docker
    # snp_dists outputs
    String? pirate_snps_dists_version = select_first([core_snp_dists.snp_dists_version,pan_snp_dists.snp_dists_version,""])
    # iqtree outputs
    String? pirate_iqtree_version = select_first([core_iqtree.version,pan_iqtree.version,""])
    # reorder matrix outputs
    File? pirate_core_snp_matrix = core_reorder_matrix.ordered_matrix
    File? pirate_iqtree_core_tree = core_reorder_matrix.tree
    File? pirate_pan_snp_matrix = pan_reorder_matrix.ordered_matrix
    File? pirate_iqtree_pan_tree = pan_reorder_matrix.tree
    # Data summary outputs
    File? pirate_summarized_data = summarize_data.summarized_data
    File? pirate_filtered_metadata = summarize_data.filtered_metadata
  }
}