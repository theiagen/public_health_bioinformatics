version 1.0

import "../../tasks/phylogenetic_inference/task_gubbins.wdl" as gubbins_task
import "../../tasks/phylogenetic_inference/task_iqtree2.wdl" as iqtree2_task
import "../../tasks/phylogenetic_inference/utilities/task_reorder_matrix.wdl" as reorder_matrix_task
import "../../tasks/phylogenetic_inference/utilities/task_snippy_core.wdl" as snippy_core_task
import "../../tasks/phylogenetic_inference/utilities/task_snp_dists.wdl" as snp_dists_task
import "../../tasks/utilities/file_handling/task_cat_files.wdl" as file_handling
import "../../tasks/phylogenetic_inference/utilities/task_shared_variants.wdl" as shared_variants_task
import "../../tasks/phylogenetic_inference/utilities/task_snp_sites.wdl" as snp_sites_task
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/data_handling/task_summarize_data.wdl" as data_summary

workflow snippy_tree_wf {
  meta {
    description: "Perform phylogenetic tree inference using iqtree (default)"
  }
  input {
    String tree_name
    Array[File] snippy_variants_outdir_tarball
    Array[String] samplenames
    File reference_genome_file
    Boolean use_gubbins = true
    Boolean core_genome = true
    Boolean call_shared_variants = true
    Array[File]? snippy_variants_qc_metrics
    
    String? data_summary_terra_project
    String? data_summary_terra_workspace
    String? data_summary_terra_table
    String? data_summary_column_names # comma delimited
    Boolean phandango_coloring = false

    # the following parameters are exposed to allow modification in snippy_streamline
    String? snippy_core_docker
    Int? snippy_core_cpu 
    Int? snippy_core_disk_size
    Int? snippy_core_memory
    File? snippy_core_bed
    
    Int? gubbins_disk_size
    Int? gubbins_memory
    Int? gubbins_cpu
    String? gubbins_docker
    
    Int? iqtree2_cpu
    Int? iqtree2_memory
    Int? iqtree2_disk_size
    String? iqtree2_opts
    String? iqtree2_docker
    Int? iqtree2_ultrafast_bootstraps
    String? iqtree2_model
    
    String? snp_dists_docker
    
    Int? snp_sites_cpu
    Int? snp_sites_disk_size
    Int? snp_sites_memory
    String? snp_sites_docker

    Boolean midpoint_root_tree = true # by default midpoint root the tree
  }
  String tree_name_updated = sub(tree_name, " ", "_")
  # snippy core creates a whole-genome multiple sequence alignment (MSA) from each alignment provided by snipy_variants
  # snippy_core does NOT create a core genome alignment- the name is misleading!
  call snippy_core_task.snippy_core {
    input:
      snippy_variants_outdir_tarball = snippy_variants_outdir_tarball,
      samplenames = samplenames,
      reference_genome_file = reference_genome_file,
      tree_name = tree_name_updated,
      docker = snippy_core_docker,
      cpu = snippy_core_cpu,
      disk_size = snippy_core_disk_size,
      memory = snippy_core_memory,
      bed_file = snippy_core_bed
  }
  # removes recombination from the MSA, if use_gubbins is set to true
  # output is a whole genome MSA without recombinant sites
  if (use_gubbins) {
    call gubbins_task.gubbins {
      input:
        alignment = snippy_core.snippy_full_alignment_clean,
        cluster_name = tree_name_updated,
        docker = gubbins_docker,
        disk_size = gubbins_disk_size,
        memory = gubbins_memory,
        cpu = gubbins_cpu
    }
  }
  # removes accessory genome sites from the MSA, creating the core genome, if core_genome is set to true
  # output will be a core genome, with or without recombinant sites removed, depending on user inputs for use_gubbins
  if (core_genome) {
    call snp_sites_task.snp_sites as snp_sites {
      input:
        # input is either the whole genome MSA, this MSA with the recombinant sites removed, 
        # or the MSA of only core sites (with or without recombinant sites as specified by use_gubbins)
        msa_fasta = select_first([gubbins.gubbins_polymorphic_fasta, snippy_core.snippy_full_alignment_clean]),
        output_name = tree_name_updated,
        output_multifasta = true,
        allow_wildcard_bases = false,
        docker = snp_sites_docker,
        output_vcf = false,
        output_phylip = false,
        output_pseudo_ref = false,
        output_monomorphic = false,
        cpu = snp_sites_cpu,
        memory = snp_sites_memory,
        disk_size = snp_sites_disk_size
    }
  }
  # creates a phylogenetic tree from the final MSA
  call iqtree2_task.iqtree2 {
    input:
      # input MSA will depend on the user-specified optional inputs for use_gubbins and core_genome
      alignment = select_first([snp_sites.snp_sites_multifasta, gubbins.gubbins_polymorphic_fasta, snippy_core.snippy_full_alignment_clean]),
      cluster_name = tree_name_updated,
      iqtree2_model = iqtree2_model,
      iqtree2_opts = iqtree2_opts,
      iqtree2_ultrafast_bootstraps = iqtree2_ultrafast_bootstraps,
      docker = iqtree2_docker,
      cpu = iqtree2_cpu,
      memory = iqtree2_memory,
      disk_size = iqtree2_disk_size
  }
  # creates a pairwise snp-distance matrix from the whole-genome MSA, with or without recombination removal.
  # whole-genome SNP matrix will always be produced regardless of whether core_genome is used 
  # because this is always valuable for interpreting strain-relatedness
  call snp_dists_task.snp_dists as wg_snp_dists {
    input:
      alignment = select_first([gubbins.gubbins_polymorphic_fasta, snippy_core.snippy_full_alignment_clean]),
      cluster_name = tree_name_updated,
      docker = snp_dists_docker
  }
  # mid-point roots the phylogenetic tree, and reorders the columns in the wgSNP matrix according to the tree tip order
  # NB the tree will remain a core genome tree is core_genome = true, and a whole-genome tree if core_genome = false
  call reorder_matrix_task.reorder_matrix as wg_reorder_matrix {
    input:
      input_tree = iqtree2.ml_tree,
      matrix = wg_snp_dists.snp_matrix,
      cluster_name = tree_name_updated + "_wg",
      midpoint_root_tree = midpoint_root_tree,
      phandango_coloring = phandango_coloring
  }
  # creates a pairwise snp-distance matrix from the core-genome MSA, if core_genome is used
  if (core_genome) {
    call snp_dists_task.snp_dists as cg_snp_dists {
      input:
        alignment = select_first([snp_sites.snp_sites_multifasta]),
        cluster_name = tree_name_updated,
        docker = snp_dists_docker
    }
    # reorders the columns in the cgSNP matrix according to the tree tip order
    # input tree is the midpoint rooted tree from the wg_reorder_matrix task, and midpoint rooting is turned off here, so the tree remains unchanged
    call reorder_matrix_task.reorder_matrix as cg_reorder_matrix {
      input:
        input_tree = wg_reorder_matrix.tree,
        matrix = cg_snp_dists.snp_matrix,
        cluster_name = tree_name_updated + "_cg",
        midpoint_root_tree = false,
        phandango_coloring = phandango_coloring
    }
  }
  # creates a data summary from comma-separated lists within single Terra data table columns
  if (defined(data_summary_column_names)) {
    call data_summary.summarize_data {
      input:
        sample_names = samplenames,
        terra_project = data_summary_terra_project,
        terra_workspace = data_summary_terra_workspace,
        terra_table = data_summary_terra_table,
        column_names = data_summary_column_names,
        output_prefix = tree_name_updated,
        phandango_coloring = phandango_coloring
    }
  }
  if (call_shared_variants) {
    call file_handling.cat_variants as concatenate_variants {
      input:
        variants_to_cat = snippy_core.snippy_variants_csv, 
        samplenames = samplenames,
        concatenated_file_name = tree_name_updated
    }
    call shared_variants_task.shared_variants {
      input:
        concatenated_variants = concatenate_variants.concatenated_variants, 
        concatenated_file_name = tree_name_updated
    }
  }
  if (defined(snippy_variants_qc_metrics)) {
    call file_handling.cat_files as concatenate_qc_metrics {
      input:
        files_to_cat = select_first([snippy_variants_qc_metrics]),
        concatenated_file_name = tree_name_updated + "_combined_qc_metrics.tsv",
        skip_extra_headers = true
    }
  }
  call versioning.version_capture {
    input:
  }
  output {
    # version capture
    String snippy_tree_version = version_capture.phb_version
    String snippy_tree_analysis_date = version_capture.date

    ### snippy core outputs
    String snippy_tree_snippy_version = snippy_core.snippy_core_version
    String snippy_tree_snippy_docker = snippy_core.snippy_core_docker_image
    File snippy_ref = snippy_core.snippy_ref
    File snippy_msa_snps_summary = snippy_core.snippy_txt

    # gubbins outputs
    String? snippy_gubbins_version = gubbins.gubbins_version
    String? snippy_gubbins_docker = gubbins.gubbins_docker
    File? snippy_gubbins_recombination_gff = gubbins.gubbins_recombination_gff
    File? snippy_gubbins_branch_stats = gubbins.gubbins_branch_stats

    ### snp_sites outputs
    String? snippy_snp_sites_version = snp_sites.snp_sites_version
    String? snippy_snp_sites_docker = snp_sites.snp_sites_docker

    ### iqtree2 outputs
    String snippy_iqtree2_version = iqtree2.iqtree2_version
    String snippy_iqtree2_docker = iqtree2.iqtree2_docker
    String snippy_iqtree2_model_used = iqtree2.iqtree2_model_used

    # snp matrix outputs
    String snippy_snp_dists_version = wg_snp_dists.snp_dists_version
    String snippy_snp_dists_docker = wg_snp_dists.snp_dists_docker
    File snippy_wg_snp_matrix = wg_reorder_matrix.ordered_matrix
    File? snippy_cg_snp_matrix = cg_reorder_matrix.ordered_matrix
    
    File snippy_final_tree = select_first([cg_reorder_matrix.tree, wg_reorder_matrix.tree]) # depending on user input for core_genome

    # data summary outputs
    File? snippy_summarized_data = summarize_data.summarized_data
    File? snippy_filtered_metadata = summarize_data.filtered_metadata

    # set final alignment from 3 possible task outputs
    File snippy_final_alignment = select_first([snp_sites.snp_sites_multifasta, gubbins.gubbins_polymorphic_fasta, snippy_core.snippy_full_alignment_clean])

    # shared snps outputs
    File? snippy_concatenated_variants = concatenate_variants.concatenated_variants
    File? snippy_shared_variants_table = shared_variants.shared_variants_table

    # combined qc metrics
    File? snippy_combined_qc_metrics = concatenate_qc_metrics.concatenated_files
  }
}
