version 1.0

import "../../utilities/task_file_handling.wdl" as file_handling
import "../../utilities/task_augur_utilities.wdl" as augur_utils
import "../../phylogenetic_inference/task_augur.wdl" as augur_tasks
import "../../phylogenetic_inference/task_snp_sites.wdl" as snp_sites_task

workflow augur {
  input {
    Array[File]+ assembly_fastas
    Array[File]+ sample_metadata_tsvs
    File reference_fasta
    File reference_genbank
    File? clades_tsv
    String build_name
    Int min_num_unambig
    Boolean use_sc2_defaults = false
  }
  call file_handling.cat_files { # concatenate all of the input fasta files together
    input:
      files_to_cat = assembly_fastas,
      concatenated_file_name = "~{build_name}_concatenated.fasta"
  }
  if (use_sc2_defaults) {
    # call get sc2 defaults
  }
  call augur_utils.filter_sequences_by_length { # remove any sequences that do not meet the quality threshold
    input:
      sequences_fasta = cat_files.concatenated_files,
      min_non_N = min_num_unambig
  }
  call augur_tasks.augur_align { # perform mafft alignment on the sequences
    input: 
      assembly_fasta = filter_sequences_by_length.filtered_fasta,
      reference_fasta = reference_fasta # make select first later
  }
  call augur_utils.tsv_join { # merge the metadata files
    input:
      input_tsvs = sample_metadata_tsvs,
      id_col = "strain",
      out_basename = "metadata-merged"
  }
  call augur_utils.fasta_to_ids { # extract list of remaining sequences (so we know which ones were dropped)
    input:
      sequences_fasta = augur_align.aligned_fasta
  }
  call snp_sites_task.snp_sites { # call variants in aligned_fasta file
    input: 
      msa_fasta = augur_align.aligned_fasta,
      output_name = build_name
  }
  call augur_tasks.augur_tree { # create a draft augur tree
    input:
      aligned_fasta = augur_align.aligned_fasta,
      build_name = build_name
  }
  call augur_tasks.augur_refine { # create a timetree (aka, refine augur tree)
    input:
      aligned_fasta = augur_align.aligned_fasta,
      draft_augur_tree = augur_tree.aligned_tree,
      metadata = tsv_join.out_tsv,
      build_name = build_name
  }
  call augur_tasks.augur_frequencies { # calculate tip frequencies
    input: 
      refined_tree = augur_refine.refined_tree,
      metadata = tsv_join.out_tsv,
      build_name = build_name
      # the following are SARS-CoV-2 specific:
      # min_date = 2020.0,
      # pivot_interval = 1,
      # pivot_interval_units = "weeks",
      # narrow_bandwidth = 0.05,
      # proportion_wide = 0.0,
  }
  call augur_tasks.augur_ancestral { # infer ancestral sequences
    input:
      refined_tree = augur_refine.refined_tree,
      aligned_fasta = augur_align.aligned_fasta,
      build_name = build_name
  }
  call augur_tasks.augur_translate { 
    input:
      refined_tree = augur_refine.refined_tree,
      ancestral_nt_muts_json = augur_ancestral.ancestral_nt_muts_json,
      reference_genbank = reference_genbank, # make into select first later
      build_name = build_name
  }
  if (defined(clades_tsv)) {
    call augur_tasks.augur_clades {
      input: 
        refined_tree = augur_refine.refined_tree,
        ancestral_nt_muts_json = augur_ancestral.ancestral_nt_muts_json,
        translated_aa_muts_json = augur_translate.translated_aa_muts_json,
        reference_fasta = reference_fasta, # make into select first
        build_name = build_name,
        clades_tsv = select_first([clades_tsv,""]) # make into select first
    }
  }
  
  
  # augur export v2 to create the auspice json file


}