version 1.0

import "../../tasks/phylogenetic_inference/task_snippy_core.wdl" as snippy_core_task
import "../../tasks/phylogenetic_inference/task_snp_sites.wdl" as snp_sites_task
import "../../tasks/phylogenetic_inference/task_iqtree2.wdl" as iqtree2_task
import "../../tasks/phylogenetic_inference/task_snp_dists.wdl" as snp_dists_task
import "../../tasks/phylogenetic_inference/task_reorder_matrix.wdl" as reorder_matrix_task
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
    File reference_genome_file
    Boolean use_gubbins = true
    Boolean core_genome = true
    String? data_summary_terra_project
    String? data_summary_terra_workspace
    String? data_summary_terra_table
    String? data_summary_column_names # comma delimited
    String? snippy_core_docker
    Int? snippy_core_cpu 
    Int? snippy_core_disk_size
    Int? snippy_core_memory
    Int? gubbins_disk_size
    Int? gubbins_memory
    Int? gubbins_cpu
    String? gubbins_docker
    Int? iqtree2_cpu
    Int? iqtree2_memory
    Int? iqtree2_disk_size
    String? iqtree2_opts
    String? iqtree2_docker
    String? iqtree2_bootstraps
    String? iqtree2_model
    String? snp_dists_docker
    Int? snp_sites_cpus
    Int? snp_sites_disk_size
    Int? snp_sites_memory
    String? snp_sites_docker
  }
  call snippy_core_task.snippy_core {
    input:
      snippy_variants_outdir_tarball = snippy_variants_outdir_tarball,
      samplenames = samplenames,
      reference_genome_file = reference_genome_file,
      tree_name = tree_name,
      docker = snippy_core_docker,
      cpu = snippy_core_cpu,
      disk_size = snippy_core_disk_size,
      memory = snippy_core_memory
  }
  if (use_gubbins) {
    call gubbins_task.gubbins {
      input:
        alignment = snippy_core.snippy_full_alignment_clean,
        cluster_name = tree_name,
        docker = gubbins_docker,
        disk_size = gubbins_disk_size,
        memory = gubbins_memory,
        cpu = gubbins_cpu
    }
  }
  if (core_genome) {
      call snp_sites_task.snp_sites as snp_sites {
        input:
          # hardcoding some of the snp-sites optional outputs to false, 
          msa_fasta = select_first([gubbins.gubbins_polymorphic_fasta,snippy_core.snippy_full_alignment_clean]),
          output_name = tree_name,
          output_multifasta = true,
          allow_wildcard_bases = false,
          docker = snp_sites_docker,
          output_vcf = false,
          output_phylip = false,
          output_pseudo_ref = false,
          output_monomorphic = false,
          cpus = snp_sites_cpus,
          memory = snp_sites_memory,
          disk_size = snp_sites_disk_size
      }
  }
  call iqtree2_task.iqtree2 {
    input:
      alignment = select_first([snp_sites.snp_sites_multifasta, gubbins.gubbins_polymorphic_fasta, snippy_core.snippy_full_alignment_clean]),
      cluster_name = tree_name,
      docker = iqtree2_docker,
      cpu = iqtree2_cpu,
      memory = iqtree2_memory,
      disk_size = iqtree2_disk_size,
      iqtree2_model = iqtree2_model,
      core_genome = core_genome,
      iqtree2_opts = iqtree2_opts,
      iqtree2_bootstraps = iqtree2_bootstraps
  }
  
  call snp_dists_task.snp_dists {
    input:
      alignment = select_first([snp_sites.snp_sites_multifasta, gubbins.gubbins_polymorphic_fasta, snippy_core.snippy_full_alignment_clean]),
      cluster_name = tree_name,
      docker = snp_dists_docker
  }
  call reorder_matrix_task.reorder_matrix {
    input:
      input_tree = iqtree2.ml_tree,
      matrix = snp_dists.snp_matrix,
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

    ### snippy core outputs ###
    String snippy_tree_snippy_version = snippy_core.snippy_version
    String snippy_tree_snippy_docker = snippy_core.snippy_docker_image
    File snippy_ref = snippy_core.snippy_ref
    File snippy_msa_snps_summary = snippy_core.snippy_txt

    # gubbins outputs
    String? snippy_gubbins_version = gubbins.version
    String? snippy_gubbins_docker = gubbins.gubbins_docker
    File? snippy_gubbins_recombination_gff = gubbins.gubbins_recombination_gff
    File? snippy_gubbins_branch_stats = gubbins.gubbins_branch_stats

    ### snp_sites outputs ###
    String? snippy_snp_sites_version = snp_sites.snp_sites_version
    String? snippy_snp_sites_docker = snp_sites.snp_sites_docker

    ### iqtree2 outputs ###
    String snippy_iqtree2_version = iqtree2.version
    String snippy_iqtree2_docker = iqtree2.iqtree2_docker
    String snippy_iqtree2_model_used = iqtree2.iqtree2_model_used

    # snpdists outputs
    String snippy_snp_dists_version = snp_dists.snp_dists_version
    String snippy_snp_dists_docker = snp_dists.snp_dists_docker

    # reorder matrix outputs
    File snippy_snp_matrix = reorder_matrix.ordered_matrix
    File snippy_final_tree = reorder_matrix.tree # this is same output tree from iqtree2, but it is midpoint rooted 

    # data summary outputs
    File? snippy_summarized_data = summarize_data.summarized_data

    # set final alignment from 3 possible task outputs
    File snippy_final_alignment = select_first([snp_sites.snp_sites_multifasta,gubbins.gubbins_polymorphic_fasta,snippy_core.snippy_full_alignment_clean])
  }
}