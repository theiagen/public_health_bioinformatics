version 1.0

import "../../tasks/phylogenetic_inference/augur/task_augur_align.wdl" as align_task
import "../../tasks/phylogenetic_inference/augur/task_augur_ancestral.wdl" as ancestral_task
import "../../tasks/phylogenetic_inference/augur/task_augur_clades.wdl" as clades_task
import "../../tasks/phylogenetic_inference/augur/task_augur_export.wdl" as export_task
import "../../tasks/phylogenetic_inference/augur/task_augur_refine.wdl" as refine_task
import "../../tasks/phylogenetic_inference/augur/task_augur_traits.wdl" as traits_task
import "../../tasks/phylogenetic_inference/augur/task_augur_translate.wdl" as translate_task
import "../../tasks/phylogenetic_inference/augur/task_augur_tree.wdl" as tree_task
import "../../tasks/phylogenetic_inference/augur/task_augur_mutation_context.wdl" as mutation_context_task
import "../../tasks/phylogenetic_inference/augur/task_extract_clade_mutations.wdl" as extract_clade_mutations_task

import "../../tasks/phylogenetic_inference/utilities/task_reorder_matrix.wdl" as reorder_matrix_task
import "../../tasks/phylogenetic_inference/utilities/task_snp_dists.wdl" as snp_dists_task

import "../../tasks/task_versioning.wdl" as versioning

import "../../tasks/utilities/data_handling/task_augur_utilities.wdl" as augur_utils
import "../../tasks/utilities/file_handling/task_cat_files.wdl" as file_handling

import "../utilities/wf_organism_parameters.wdl" as set_organism_defaults

workflow augur {
  input {
    # primary inputs
    Array[File]+ assembly_fastas # use the HA or NA segment files for flu
    Array[File]? sample_metadata_tsvs # created with Augur_Prep
    String build_name
    File? alignment_fasta # if alignment is provided, then skip alignment step

    # organism-specific inputs
    String organism = "sars-cov-2" # compatible with organism_parameters inputs, or manual Augur parameters below
    String flu_segment = "HA" # options: HA or NA
    String? flu_subtype # options: "Victoria" "Yamagata" "H3N2" "H1N1" "H5N1"

    # augur parameters
    File? reference_fasta
    Boolean remove_reference = false # by default, do not remove the reference
    File? reference_genbank
    Int? min_num_unambig
    File? clades_tsv
    File? lat_longs_tsv
    File? auspice_config
    Int? pivot_interval
    Float? min_date
    Float? narrow_bandwidth
    Float? proportion_wide
    String? augur_trait_columns # comma-separated list of columns to use for traits
    Boolean generate_clades_tsv = false # generate clades tsv file from clade_membership header
    String augur_id_column = "strain" # column in metadata tsv that contains the sequence names/IDs

    # phylogenetic tree = true # by default, midpoint root the tree
    Boolean midpoint_root_tree = true
    String? outgroup_root
  }
  String build_name_updated = sub(build_name, " ", "_")

  # capture the version
  call versioning.version_capture { 
    input:
  }
  # skip alignment if alignment_fasta is inputted
  if (! defined(alignment_fasta)) {
    Boolean call_alignment = true
  }

  # set organism parameters for default organisms, passthrough for others
  call set_organism_defaults.organism_parameters {
    input:
      organism = organism,
      reference_genbank = reference_genbank,
      reference_genome = reference_fasta,
      min_num_unambig = min_num_unambig,
      flu_segment = flu_segment,
      flu_subtype = flu_subtype,
      clades_tsv = clades_tsv,
      lat_longs_tsv = lat_longs_tsv,
      auspice_config = auspice_config,
      pivot_interval = pivot_interval,
      min_date = min_date,
      narrow_bandwidth = narrow_bandwidth,
      proportion_wide = proportion_wide
  }
  # skip clade extraction if augur_clade_columns is not defined
  if (defined(clades_tsv) || (defined(organism_parameters.augur_clades_tsv) && (basename(organism_parameters.augur_clades_tsv) != "minimal-clades.tsv"))) {
    Boolean call_clades = true
  } 
  if (defined(sample_metadata_tsvs)) {
    # merge the metadata files
    call augur_utils.tsv_join { 
      input:
        input_tsvs = select_first([sample_metadata_tsvs]),
        id_col = augur_id_column,
        out_basename = "metadata-merged"
    }
  }
  if (defined(call_alignment)) {
    # concatenate all of the input fasta files together
    call file_handling.cat_files {       
      input:
        files_to_cat = assembly_fastas,
        concatenated_file_name = "~{build_name_updated}_concatenated.fasta"
    }
  }

  # Alignment preparation
  # remove any sequences that do not meet the quality threshold
  # perform prior to alignment to increase throughput
  call augur_utils.filter_sequences_by_length { 
    input:
      sequences_fasta = select_first([cat_files.concatenated_files, alignment_fasta]),
      min_non_N = select_first([min_num_unambig, organism_parameters.augur_min_num_unambig]),
  }
  if (defined(call_alignment)) {
    # perform mafft alignment on the sequences
    call align_task.augur_align { 
      input:
        assembly_fasta = filter_sequences_by_length.filtered_fasta,
        reference_fasta = select_first([reference_fasta, organism_parameters.reference]),
        remove_reference = remove_reference
    }
  }
  call augur_utils.fasta_to_ids { 
    # extract list of remaining sequences (so we know which ones were dropped)
    input:
      sequences_fasta = select_first([augur_align.aligned_fasta, filter_sequences_by_length.filtered_fasta])
  }

  # Phylogenetic tree reconstruction
  # create a phylogenetic tree
  call tree_task.augur_tree { 
    input:
      aligned_fasta = select_first([augur_align.aligned_fasta, filter_sequences_by_length.filtered_fasta]),
      build_name = build_name_updated
  }
  # create a snp matrix from the alignment
  call snp_dists_task.snp_dists { 
    input:
      cluster_name = build_name_updated,
      alignment = select_first([augur_align.aligned_fasta, filter_sequences_by_length.filtered_fasta])
  }
  # reorder snp matrix to match distance tree 
  call reorder_matrix_task.reorder_matrix { 
    input:
      input_tree = augur_tree.tree,
      matrix = snp_dists.snp_matrix,
      cluster_name = build_name_updated,
      midpoint_root_tree = midpoint_root_tree,
      outgroup_root = outgroup_root
  }

  # refine and potentially create a time-calibrated phylogenetic tree
  call refine_task.augur_refine { 
    input:
      aligned_fasta = select_first([augur_align.aligned_fasta, filter_sequences_by_length.filtered_fasta]),
      draft_augur_tree = reorder_matrix.tree,
      metadata = tsv_join.out_tsv,
      build_name = build_name_updated,
      build_time_tree = tsv_join.has_time
  }

  if (defined(tsv_join.out_tsv)) {
    # infer ancestral sequences
    call ancestral_task.augur_ancestral { 
      input:
        refined_tree = augur_refine.refined_tree,
        aligned_fasta = select_first([augur_align.aligned_fasta, filter_sequences_by_length.filtered_fasta]),
        build_name = build_name_updated
    }
    # translate gene regions from nucleotides to amino acids
    call translate_task.augur_translate { 
      input:
        refined_tree = augur_refine.refined_tree,
        ancestral_nt_muts_json = augur_ancestral.ancestral_nt_muts_json,
        reference_genbank = select_first([reference_genbank, organism_parameters.reference_gbk]),
        build_name = build_name_updated
    }
    if (organism_parameters.standardized_organism == "MPXV") {
      call mutation_context_task.mutation_context { # add mutation context to the tree
        input:
          refined_tree = augur_refine.refined_tree,
          ancestral_nt_muts_json = augur_ancestral.ancestral_nt_muts_json,
          build_name = build_name_updated
      }
    }
    call traits_task.augur_traits {
      input:
        refined_tree = augur_refine.refined_tree,
        metadata = tsv_join.out_tsv,
        columns = select_first([augur_trait_columns, "pango_lineage,clade_membership"]), # default to these columns if none are specified
        build_name = build_name_updated
    }
    if (defined(call_clades)) {
      if (generate_clades_tsv) { 
        call extract_clade_mutations_task.extract_clade_mutations {
          input:
            metadata_tsv = select_first([tsv_join.out_tsv]),
            clade_columns = "clade_membership",
            tip_column = augur_id_column,
            tree = augur_refine.refined_tree,
            nt_mutations = augur_ancestral.ancestral_nt_muts_json,
            aa_mutations = augur_translate.translated_aa_muts_json
        }
      }
      # assign clades to nodes based on amino-acid or nucleotide signatures
      call clades_task.augur_clades { 
        input:
          refined_tree = augur_refine.refined_tree,
          ancestral_nt_muts_json = augur_ancestral.ancestral_nt_muts_json,
          translated_aa_muts_json = augur_translate.translated_aa_muts_json,
          build_name = build_name_updated,
          clades_tsv = select_first([extract_clade_mutations.clades_tsv, clades_tsv, organism_parameters.augur_clades_tsv])
      }
    }
  }
  # export json files suitable for auspice visualization
  call export_task.augur_export { 
    input:
      tree = select_first([augur_refine.refined_tree, augur_tree.tree]),
      metadata = select_first([tsv_join.out_tsv, "gs://theiagen-public-resources-rp/empty_files/empty.fasta"]),
      node_data_jsons = select_all([
                          augur_refine.branch_lengths,
                          augur_ancestral.ancestral_nt_muts_json,
                          augur_translate.translated_aa_muts_json,
                          augur_clades.clade_assignments_json,
                          augur_traits.traits_assignments_json,
                          mutation_context.mutation_context_json]),
      build_name = build_name_updated,
      lat_longs_tsv = select_first([lat_longs_tsv, organism_parameters.augur_lat_longs_tsv]),
      auspice_config = select_first([auspice_config, organism_parameters.augur_auspice_config])
  }

  # determine what the refined tree represents
  if (select_first([tsv_join.has_time, false])) {
    File time_tree_path = select_first([augur_refine.refined_tree])
  }
  if (! select_first([tsv_join.has_time, false])) {
    File phylogenetic_tree_path = select_first([augur_refine.refined_tree])
  }
  output {
    # version capture
    String augur_phb_version = version_capture.phb_version
    String augur_phb_analysis_date = version_capture.date
    String augur_version = augur_tree.augur_version

    # augur outputs
    String? augur_mafft_version = augur_align.mafft_version
    File? auspice_input_json = augur_export.auspice_json
    File? time_tree = time_tree_path
    File phylogenetic_tree = select_first([phylogenetic_tree_path, augur_tree.tree])
    String augur_iqtree_model_used = augur_tree.iqtree_model_used
    String augur_iqtree_version = augur_tree.iqtree_version
    String augur_fasttree_version = augur_tree.fasttree_version
    String augur_raxml_version = augur_tree.raxml_version
    File alignment_fasta = select_first([augur_align.aligned_fasta, alignment_fasta])
    File combined_assemblies = filter_sequences_by_length.filtered_fasta
    File? metadata_merged = tsv_join.out_tsv
    File? traits_json = augur_traits.traits_assignments_json

    # clade assignments
    File? clade_mutations = extract_clade_mutations.clades_tsv

    # list of samples that were kept and met the length filters    
    File keep_list = fasta_to_ids.ids_txt
  
    # snp matrix output
    File snp_matrix = reorder_matrix.ordered_matrix
  }
}