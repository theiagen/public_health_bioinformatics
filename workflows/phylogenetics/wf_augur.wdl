version 1.0

import "../../tasks/utilities/task_file_handling.wdl" as file_handling
import "../../tasks/utilities/task_augur_utilities.wdl" as augur_utils

import "../../tasks/phylogenetic_inference/augur/task_augur_align.wdl" as align_task
import "../../tasks/phylogenetic_inference/augur/task_augur_ancestral.wdl" as ancestral_task
import "../../tasks/phylogenetic_inference/augur/task_augur_clades.wdl" as clades_task
import "../../tasks/phylogenetic_inference/augur/task_augur_export.wdl" as export_task
import "../../tasks/phylogenetic_inference/augur/task_augur_frequencies.wdl" as frequencies_task
import "../../tasks/phylogenetic_inference/augur/task_augur_refine.wdl" as refine_task
import "../../tasks/phylogenetic_inference/augur/task_augur_translate.wdl" as translate_task
import "../../tasks/phylogenetic_inference/augur/task_augur_tree.wdl" as tree_task

import "../../tasks/phylogenetic_inference/task_snp_sites.wdl" as snp_sites_task
import "../../tasks/phylogenetic_inference/task_snp_dists.wdl" as snp_dists_task
import "../../tasks/phylogenetic_inference/task_reorder_matrix.wdl" as reorder_matrix_task

import "../../tasks/task_versioning.wdl" as versioning

workflow augur {
  input {
    Array[File]+ assembly_fastas # use ha or na segments for flu
    Array[File]+ sample_metadata_tsvs
    String build_name
    File? reference_fasta
    File? reference_genbank
    Int? min_num_unambig
    String organism = "sars-cov-2" # options: sars-cov-2 or flu
    String flu_segment = "HA" # options: HA or NA
    String? flu_subtype # options: "Victoria" "Yamagata" "H3N2" "H1N1"

    # these following inputs should be optional, but I'm worried the select_first will make them "not optional" 
    # this is especially the case if organism is not given
    File? clades_tsv
    File? lat_longs_tsv 
    File? auspice_config
    Float? min_frequency_date

    Boolean distance_tree_only = false
  }
  call file_handling.cat_files { # concatenate all of the input fasta files together
    input:
      files_to_cat = assembly_fastas,
      concatenated_file_name = "~{build_name}_concatenated.fasta"
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
  call augur_utils.filter_sequences_by_length { # remove any sequences that do not meet the quality threshold
    input:
      sequences_fasta = cat_files.concatenated_files,
      min_non_N = select_first([sc2_defaults.min_num_unambig, flu_defaults.min_num_unambig, min_num_unambig])
  }
  call align_task.augur_align { # perform mafft alignment on the sequences
    input: 
      assembly_fasta = filter_sequences_by_length.filtered_fasta,
      reference_fasta = select_first([reference_fasta, sc2_defaults.reference_fasta, flu_defaults.reference_fasta])
  }
  call augur_utils.tsv_join { # merge the metadata files
    input:
      input_tsvs = sample_metadata_tsvs,
      id_col = "strain",
      out_basename = "metadata-merged"
  }

  ## keep the following two tasks???
  call augur_utils.fasta_to_ids { # extract list of remaining sequences (so we know which ones were dropped)
    input:
      sequences_fasta = augur_align.aligned_fasta
  }
  call snp_sites_task.snp_sites { # call variants in aligned_fasta file
    input: 
      msa_fasta = augur_align.aligned_fasta,
      output_name = build_name
  }
  ## keep the preceeding two tasks???

  call tree_task.augur_tree { # create a "draft" augur tree
    input:
      aligned_fasta = augur_align.aligned_fasta,
      build_name = build_name
  }
  if (! distance_tree_only) { # by default, continue
    call refine_task.augur_refine { # create a timetree (aka, refine augur tree)
      input:
        aligned_fasta = augur_align.aligned_fasta,
        draft_augur_tree = augur_tree.aligned_tree,
        metadata = tsv_join.out_tsv,
        build_name = build_name
    }
    # call frequencies_task.augur_frequencies { # calculate tip frequencies
    #   input: 
    #     refined_tree = augur_refine.refined_tree,
    #     metadata = tsv_join.out_tsv,
    #     build_name = build_name,
    #     min_date = select_first([sc2_defaults.min_date, flu_defaults.min_date, min_frequency_date]), 
    #     pivot_interval = select_first([sc2_defaults.pivot_interval, flu_defaults.pivot_interval, 3]),
    #     pivot_interval_units = select_first([sc2_defaults.pivot_interval_units, "months"]),
    #     narrow_bandwidth = select_first([sc2_defaults.narrow_bandwidth, flu_defaults.narrow_bandwidth, 0.08333333333333333]),
    #     proportion_wide = select_first([sc2_defaults.proportion_wide, flu_defaults.proportion_wide, 0.2])
    # }
    call ancestral_task.augur_ancestral { # infer ancestral sequences
      input:
        refined_tree = augur_refine.refined_tree,
        aligned_fasta = augur_align.aligned_fasta,
        build_name = build_name
    }
    call translate_task.augur_translate { # translate gene regions from nucleotides to amino acids
      input:
        refined_tree = augur_refine.refined_tree,
        ancestral_nt_muts_json = augur_ancestral.ancestral_nt_muts_json,
        reference_genbank = select_first([reference_genbank, sc2_defaults.reference_genbank, flu_defaults.reference_genbank]),
        build_name = build_name
    }
    if (flu_segment == "ha") { # only have clade information for ha segments
      call clades_task.augur_clades { # assign clades to nodes based on amino-acid or nucleotide signatures
        input: 
          refined_tree = augur_refine.refined_tree,
          ancestral_nt_muts_json = augur_ancestral.ancestral_nt_muts_json,
          translated_aa_muts_json = augur_translate.translated_aa_muts_json,
          reference_fasta = select_first([reference_fasta, sc2_defaults.reference_fasta, flu_defaults.reference_fasta]),
          build_name = build_name,
          clades_tsv = select_first([clades_tsv, sc2_defaults.clades_tsv, flu_defaults.clades_tsv])
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
                            augur_clades.clade_assignments_json]),
        build_name = build_name,
        lat_longs_tsv = select_first([lat_longs_tsv, sc2_defaults.lat_longs_tsv, flu_defaults.lat_longs_tsv]),
        auspice_config = select_first([auspice_config, sc2_defaults.auspice_config, flu_defaults.auspice_config])
    }
  }
  call snp_dists_task.snp_dists { # create a snp matrix from the alignment
    input:
      cluster_name = build_name,
      alignment = augur_align.aligned_fasta
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
    String augur_version = augur_align.augur_version

    # augur outputs
    File? auspice_input_json = augur_export.auspice_json
    File? time_tree = augur_refine.refined_tree
    File distance_tree = augur_tree.aligned_tree
    File aligned_fastas = augur_align.aligned_fasta
    File combined_assemblies = filter_sequences_by_length.filtered_fasta
    File metadata_merged = tsv_join.out_tsv
    
    # not sure if wanting to keep the tasks that make these
    File keep_list = fasta_to_ids.ids_txt
    File? unmasked_snps = snp_sites.snp_sites_vcf
  
    # snp matrix output
    File snp_matrix = reorder_matrix.ordered_matrix
  }
}