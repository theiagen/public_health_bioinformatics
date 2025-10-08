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

import "../../tasks/phylogenetic_inference/utilities/task_reorder_matrix.wdl" as reorder_matrix_task
import "../../tasks/phylogenetic_inference/utilities/task_snp_dists.wdl" as snp_dists_task

import "../../tasks/task_versioning.wdl" as versioning

import "../../tasks/utilities/data_handling/task_augur_utilities.wdl" as augur_utils
import "../../tasks/utilities/file_handling/task_cat_files.wdl" as file_handling

import "../utilities/wf_organism_parameters.wdl" as set_organism_defaults

import "../../tasks/phylogenetic_inference/utilities/task_root_phylo.wdl" as task_root_phylo

workflow augur {
  input {
    Array[File]+ assembly_fastas # use the HA or NA segment files for flu
    Array[File]? sample_metadata_tsvs # created with Augur_Prep
    String build_name
    String build_name_updated = sub(build_name, " ", "_")
    File? reference_fasta
    Boolean remove_reference = false # by default, do not remove the reference
    File? reference_genbank
    Int? min_num_unambig
    String organism = "sars-cov-2" # options: sars-cov-2, flu, mpxv, "rsv-a" or "rsv-b"
    String flu_segment = "HA" # options: HA or NA
    String? flu_subtype # options: "Victoria" "Yamagata" "H3N2" "H1N1" "H5N1"
    Boolean skip_alignment = false # by default, do not skip alignment
    File? alignment_fasta # if alignment is skipped, provide an alignment

    File? clades_tsv
    Boolean extract_clade_mutations = false # generate clades_tsv on the fly
    Boolean run_traits = false # by default, do not run traits
    String? augur_trait_columns # comma-separated list of columns to use for traits
    # these are very minimal files that hopefully will prevent workflow failure but will not provide any useful information
    File? lat_longs_tsv
    File? auspice_config

    Boolean distance_tree_only = false # by default, do not skip making a time tree

    Boolean midpoint_root_tree = true # by default, midpoint root the tree
    String? outgroup_root

    Int? pivot_interval
    Float? min_date
    Float? narrow_bandwidth
    Float? proportion_wide
  }
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
  if (defined(sample_metadata_tsvs)) {
    call augur_utils.tsv_join { # merge the metadata files
      input:
        input_tsvs = select_first([sample_metadata_tsvs]),
        id_col = "strain",
        out_basename = "metadata-merged"
    }
  }
  if (! skip_alignment) { # by default, continue
    call file_handling.cat_files { # concatenate all of the input fasta files together
      input:
        files_to_cat = assembly_fastas,
        concatenated_file_name = "~{build_name_updated}_concatenated.fasta"
    }
  }
  call augur_utils.filter_sequences_by_length { # remove any sequences that do not meet the quality threshold
    input:
      sequences_fasta = select_first([cat_files.concatenated_files, alignment_fasta]),
      min_non_N = select_first([min_num_unambig, organism_parameters.augur_min_num_unambig]),
  }
  if (! skip_alignment) { # by default, continue
    call align_task.augur_align { # perform mafft alignment on the sequences
      input:
        assembly_fasta = filter_sequences_by_length.filtered_fasta,
        reference_fasta = select_first([reference_fasta, organism_parameters.reference]),
        remove_reference = remove_reference
    }
  }
  call augur_utils.fasta_to_ids { # extract list of remaining sequences (so we know which ones were dropped)
    input:
      sequences_fasta = select_first([augur_align.aligned_fasta, filter_sequences_by_length.filtered_fasta])
  }
  call tree_task.augur_tree { # create a "draft" (or distance) augur tree
    input:
      aligned_fasta = select_first([augur_align.aligned_fasta, filter_sequences_by_length.filtered_fasta]),
      build_name = build_name_updated
  }
  if (! distance_tree_only && defined(tsv_join.out_tsv)) { # by default, continue
    call refine_task.augur_refine { # create a timetree (aka, refine augur tree)
      input:
        aligned_fasta = select_first([augur_align.aligned_fasta, filter_sequences_by_length.filtered_fasta]),
        draft_augur_tree = augur_tree.aligned_tree,
        metadata = tsv_join.out_tsv,
        build_name = build_name_updated
    }
    call ancestral_task.augur_ancestral { # infer ancestral sequences
      input:
        refined_tree = augur_refine.refined_tree,
        aligned_fasta = select_first([augur_align.aligned_fasta, filter_sequences_by_length.filtered_fasta]),
        build_name = build_name_updated
    }
    call translate_task.augur_translate { # translate gene regions from nucleotides to amino acids
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
    if (flu_segment != "NA") { # we have clade information for all "standard" species except for NA flu segments (SC2 defaults should be selected first)
      if (run_traits && defined(tsv_join.out_tsv)) { # by default do not run traits and clades will be assigned based on the clades_tsv
        call traits_task.augur_traits {
          input:
            refined_tree = augur_refine.refined_tree,
            metadata = tsv_join.out_tsv,
            columns = select_first([augur_trait_columns, "pango_lineage,clade_membership"]), # default to these columns if none are specified
            build_name = build_name_updated
        }
      }
      if (! run_traits) {
        if (defined(clades_tsv) || (defined(organism_parameters.augur_clades_tsv) && (basename(organism_parameters.augur_clades_tsv) != "minimal-clades.tsv"))) { # one must be present and not the empty "minimal-clades.tsv" file
          call clades_task.augur_clades { # assign clades to nodes based on amino-acid or nucleotide signatures
            input:
              refined_tree = augur_refine.refined_tree,
              ancestral_nt_muts_json = augur_ancestral.ancestral_nt_muts_json,
              translated_aa_muts_json = augur_translate.translated_aa_muts_json,
              build_name = build_name_updated,
              clades_tsv = select_first([clades_tsv, organism_parameters.augur_clades_tsv])
          }
        }
      }
    }
    call export_task.augur_export { # export json files suitable for auspice visualization
      input:
        refined_tree = augur_refine.refined_tree,
        metadata = tsv_join.out_tsv,
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
  }
  call snp_dists_task.snp_dists { # create a snp matrix from the alignment
    input:
      cluster_name = build_name_updated,
      alignment = select_first([augur_align.aligned_fasta,filter_sequences_by_length.filtered_fasta])
  }
  call reorder_matrix_task.reorder_matrix { # reorder snp matrix to match distance tree 
    input:
      input_tree = augur_tree.aligned_tree,
      matrix = snp_dists.snp_matrix,
      cluster_name = build_name_updated,
      midpoint_root_tree = midpoint_root_tree
  }
  call versioning.version_capture { # capture the version
    input:
  }
  output {
    # version capture
    String augur_phb_version = version_capture.phb_version
    String augur_phb_analysis_date = version_capture.date
    String augur_version = augur_tree.augur_version

    # augur outputs
    String? augur_mafft_version = augur_align.mafft_version
    File? auspice_input_json = augur_export.auspice_json
    File? time_tree = augur_refine.refined_tree
    File distance_tree = augur_tree.aligned_tree
    String augur_iqtree_model_used = augur_tree.iqtree_model_used
    String augur_iqtree_version = augur_tree.iqtree_version
    String augur_fasttree_version = augur_tree.fasttree_version
    String augur_raxml_version = augur_tree.raxml_version
    File aligned_fastas = select_first([augur_align.aligned_fasta, alignment_fasta])
    File combined_assemblies = filter_sequences_by_length.filtered_fasta
    File? metadata_merged = tsv_join.out_tsv
    File? traits_json = augur_traits.traits_assignments_json

    # list of samples that were kept and met the length filters    
    File keep_list = fasta_to_ids.ids_txt
  
    # snp matrix output
    File snp_matrix = reorder_matrix.ordered_matrix
  }
}