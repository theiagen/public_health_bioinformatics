version 1.0

import "../../tasks/utilities/task_file_handling.wdl" as file_handling
import "../../tasks/utilities/task_augur_utilities.wdl" as augur_utils

import "../../tasks/phylogenetic_inference/augur/task_augur_align.wdl" as align_task
import "../../tasks/phylogenetic_inference/augur/task_augur_ancestral.wdl" as ancestral_task
import "../../tasks/phylogenetic_inference/augur/task_augur_traits.wdl" as traits_task
import "../../tasks/phylogenetic_inference/augur/task_augur_clades.wdl" as clades_task
import "../../tasks/phylogenetic_inference/augur/task_augur_export.wdl" as export_task
import "../../tasks/phylogenetic_inference/augur/task_augur_refine.wdl" as refine_task
import "../../tasks/phylogenetic_inference/augur/task_augur_translate.wdl" as translate_task
import "../../tasks/phylogenetic_inference/augur/task_augur_tree.wdl" as tree_task

import "../../tasks/phylogenetic_inference/task_snp_dists.wdl" as snp_dists_task
import "../../tasks/phylogenetic_inference/task_reorder_matrix.wdl" as reorder_matrix_task

import "../../tasks/task_versioning.wdl" as versioning

workflow augur {
  input {
    Array[File]+ assembly_fastas # use the HA or NA segment files for flu
    Array[File]+ sample_metadata_tsvs # created with Augur_Prep
    String build_name
    File? reference_fasta
    File? reference_genbank
    Int? min_num_unambig
    String organism = "sars-cov-2" # options: sars-cov-2 or flu or mpxv
    String flu_segment = "HA" # options: HA or NA
    String? flu_subtype # options: "Victoria" "Yamagata" "H3N2" "H1N1"
    Boolean skip_alignment = false # by default, do not skip alignment
    File? alignment_fasta # if alignment is skipped, provide an alignment

    File? clades_tsv
    Boolean run_traits = false # by default, do not run traits
    String? augur_trait_columns # comma-separated list of columns to use for traits
    # these are very minimal files that hopefully will prevent workflow failure but will not provide any useful information
    File lat_longs_tsv = "gs://theiagen-public-files-rp/terra/augur-defaults/minimal-lat-longs.tsv"
    File auspice_config = "gs://theiagen-public-files-rp/terra/augur-defaults/minimal-auspice-config.json"
    File genes = "gs://theiagen-public-files-rp/terra/augur-defaults/minimal-genes.tsv"
    File colors = "gs://theiagen-public-files-rp/terra/augur-defaults/minimal-colors.tsv"

    Boolean distance_tree_only = false # by default, do not skip making a time tree
  }
  if (organism == "sars-cov-2") {
    call augur_utils.set_sc2_defaults as sc2_defaults { # establish default parameters for sars-cov-2
      input:
    }
  }
  if (organism == "flu") {
    call augur_utils.set_flu_defaults as flu_defaults { # establish default parameters for flu
      input:
        flu_segment = flu_segment,
        flu_subtype = flu_subtype
    }
  }
  if (organism == "MPXV" || organism == "mpxv" || organism == "monkeypox") {
    call augur_utils.set_mpxv_defaults as mpxv_defaults { # establish default parameters for mpxv
      input:
    }
  }
  call augur_utils.tsv_join { # merge the metadata files
    input:
      input_tsvs = sample_metadata_tsvs,
      id_col = "strain",
      out_basename = "metadata-merged"
  }
  if (! skip_alignment) { # by default, continue
    call file_handling.cat_files { # concatenate all of the input fasta files together
      input:
        files_to_cat = assembly_fastas,
        concatenated_file_name = "~{build_name}_concatenated.fasta"
    }
  }
  call augur_utils.filter_sequences_by_length { # remove any sequences that do not meet the quality threshold
    input:
      sequences_fasta = select_first([cat_files.concatenated_files, alignment_fasta]),
      min_non_N = select_first([min_num_unambig, sc2_defaults.min_num_unambig, flu_defaults.min_num_unambig, mpxv_defaults.min_num_unambig]),
  }
  if (! skip_alignment) { # by default, continue
    call align_task.augur_align { # perform mafft alignment on the sequences
      input:
        assembly_fasta = filter_sequences_by_length.filtered_fasta,
        reference_fasta = select_first([reference_fasta, sc2_defaults.reference_fasta, flu_defaults.reference_fasta, mpxv_defaults.reference_fasta]),
    }
  }
  call augur_utils.fasta_to_ids { # extract list of remaining sequences (so we know which ones were dropped)
    input:
      sequences_fasta = select_first([augur_align.aligned_fasta, filter_sequences_by_length.filtered_fasta])
  }
  call tree_task.augur_tree { # create a "draft" (or distance) augur tree
    input:
      aligned_fasta = select_first([augur_align.aligned_fasta, filter_sequences_by_length.filtered_fasta]),
      build_name = build_name
  }
  if (! distance_tree_only) { # by default, continue
    call refine_task.augur_refine { # create a timetree (aka, refine augur tree)
      input:
        aligned_fasta = select_first([augur_align.aligned_fasta, filter_sequences_by_length.filtered_fasta]),
        draft_augur_tree = augur_tree.aligned_tree,
        metadata = tsv_join.out_tsv,
        build_name = build_name
    }
    call ancestral_task.augur_ancestral { # infer ancestral sequences
      input:
        refined_tree = augur_refine.refined_tree,
        aligned_fasta = select_first([augur_align.aligned_fasta, filter_sequences_by_length.filtered_fasta]),
        build_name = build_name
    }
    call translate_task.augur_translate { # translate gene regions from nucleotides to amino acids
      input:
        refined_tree = augur_refine.refined_tree,
        ancestral_nt_muts_json = augur_ancestral.ancestral_nt_muts_json,
        reference_genbank = select_first([reference_genbank, sc2_defaults.reference_genbank, flu_defaults.reference_genbank, mpxv_defaults.reference_genbank]),
        genes = select_first([mpxv_defaults.genes, genes]),
        build_name = build_name
    }
    if (flu_segment == "HA") { # we only have clade information for HA segments (but SC2 defaults will be selected first)
      if (run_traits) { # by default do not run traits and clades will be assigned based on the clades_tsv
        call traits_task.augur_traits {
          input:
            refined_tree = augur_refine.refined_tree,
            metadata = tsv_join.out_tsv,
            columns = select_first([augur_trait_columns, "pango_lineage,clade_membership"]), # default to these columns if none are specified
            build_name = build_name
        }
      }
      if (! run_traits) {
        if (defined(clades_tsv) || defined(sc2_defaults.clades_tsv) || defined(flu_defaults.clades_tsv) || defined(mpxv_defaults.clades_tsv)) { # one of these must be present
          call clades_task.augur_clades { # assign clades to nodes based on amino-acid or nucleotide signatures
            input:
              refined_tree = augur_refine.refined_tree,
              ancestral_nt_muts_json = augur_ancestral.ancestral_nt_muts_json,
              translated_aa_muts_json = augur_translate.translated_aa_muts_json,
              reference_fasta = select_first([reference_fasta, sc2_defaults.reference_fasta, flu_defaults.reference_fasta, mpxv_defaults.reference_fasta]),
              build_name = build_name,
              clades_tsv = select_first([clades_tsv, sc2_defaults.clades_tsv, flu_defaults.clades_tsv, mpxv_defaults.clades_tsv])
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
                            augur_traits.traits_assignments_json]),
        build_name = build_name,
        colors_tsv = select_first([mpxv_defaults.colors, colors]),
        lat_longs_tsv = select_first([sc2_defaults.lat_longs_tsv, flu_defaults.lat_longs_tsv, mpxv_defaults.lat_longs_tsv, lat_longs_tsv]),
        auspice_config = select_first([sc2_defaults.auspice_config, flu_defaults.auspice_config, mpxv_defaults.auspice_config, auspice_config])
    }
  }
  call snp_dists_task.snp_dists { # create a snp matrix from the alignment
    input:
      cluster_name = build_name,
      alignment = select_first([augur_align.aligned_fasta,filter_sequences_by_length.filtered_fasta])
  }
  call reorder_matrix_task.reorder_matrix { # reorder snp matrix to match distance tree 
    input:
      input_tree = augur_tree.aligned_tree,
      matrix = snp_dists.snp_matrix,
      cluster_name = build_name
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
    File? auspice_input_json = augur_export.auspice_json
    File? time_tree = augur_refine.refined_tree
    File distance_tree = augur_tree.aligned_tree
    File aligned_fastas = select_first([augur_align.aligned_fasta, alignment_fasta])
    File combined_assemblies = filter_sequences_by_length.filtered_fasta
    File metadata_merged = tsv_join.out_tsv
    File? traits_json = augur_traits.traits_assignments_json

    # list of samples that were kept and met the length filters    
    File keep_list = fasta_to_ids.ids_txt
  
    # snp matrix output
    File snp_matrix = reorder_matrix.ordered_matrix
  }
}