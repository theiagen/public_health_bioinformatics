version 1.0

import "../tasks/tasks_nextstrain.wdl" as nextstrain
import "../tasks/tasks_utils.wdl" as utils
import "../tasks/task_phylo.wdl" as phylo
import "../tasks/task_versioning.wdl" as versioning

workflow theiacov_distance_tree {
  meta {
    description: "Workflow for SC2 cluster investigations. TheiaCoV_Augur_DistanceTree is will generate a ML distance tree using select tasks incorporated in the ThieaCoV_Augur_Run workflow; output from the modified sarscov2_nextstrain workflow will also be used to infer SNP distances. The ML distance tree output can be visualized using the Auspice web application https://auspice.us/"
    author: "Kevin G Libuit"
    email:  "kevin.libuit@theiagen.com"
  }
  input {
    Array[File]+ assembly_fastas
    Array[File]+ sample_metadata_tsvs
    String build_name
    File? builds_yaml
    File? ref_fasta
    Int min_unambig_genome = 27000
  }
  parameter_meta {
    assembly_fastas: {
    description: "Set of assembled genomes to align and build trees. These must represent a single chromosome/segment of a genome only. Fastas may be one-sequence-per-individual or a concatenated multi-fasta (unaligned) or a mixture of the two. They may be compressed (gz, bz2, zst, lz4), uncompressed, or a mixture.",
    patterns: ["*.fasta", "*.fa", "*.fasta.gz", "*.fasta.zst"]
    }
    sample_metadata_tsvs: {
    description: "Tab-separated metadata file that contain binning variables and values. Must contain all samples: output will be filtered to the IDs present in this file.",
    patterns: ["*.txt", "*.tsv"]
    }
    ref_fasta: {
    description: "A reference assembly (not included in assembly_fastas) to align assembly_fastas against. Typically from NCBI RefSeq or similar.",
    patterns: ["*.fasta", "*.fa"]
    }
    min_unambig_genome: {
    description: "Minimum number of called bases in genome to pass prefilter."
    }
    ancestral_traits_to_infer: {
    description: "A list of metadata traits to use for ancestral node inference (see https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/traits.html). Multiple traits may be specified; must correspond exactly to column headers in metadata file. Omitting these values will skip ancestral trait inference, and ancestral nodes will not have estimated values for metadata."
    }
    clades_tsv: {
    description: "A TSV file containing clade mutation positions in four columns: [clade  gene    site    alt]; see: https://nextstrain.org/docs/tutorials/defining-clades",
    patterns: ["*.tsv", "*.txt"]
    }
  }
  call nextstrain.nextstrain_ncov_defaults
  #### mafft_and_snp
  call utils.zcat {
    input:
      infiles = assembly_fastas,
      output_name = "all_samples_combined_assembly.fasta"
    }
  call nextstrain.nextstrain_deduplicate_sequences as dedup_seqs {
    input:
      sequences_fasta = zcat.combined
    }
  call utils.filter_sequences_by_length {
    input:
      sequences_fasta = dedup_seqs.sequences_deduplicated_fasta,
      min_non_N = min_unambig_genome
    }
  call nextstrain.mafft_one_chr_chunked as mafft {
    input:
      sequences = filter_sequences_by_length.filtered_fasta,
      ref_fasta = select_first([ref_fasta, nextstrain_ncov_defaults.reference_fasta]),
      basename = "all_samples_aligned.fasta"
    }
  #### merge metadata, compute derived cols
  if(length(sample_metadata_tsvs)>1) {
    call utils.tsv_join {
      input:
        input_tsvs = sample_metadata_tsvs,
        id_col = 'strain',
        out_basename = "metadata-merged"
    }
  }
  call nextstrain.derived_cols {
    input:
      metadata_tsv = select_first(flatten([[tsv_join.out_tsv], sample_metadata_tsvs]))
  }
  ## Subsample if builds.yaml file provided
  if(defined(builds_yaml)) {
    call nextstrain.nextstrain_build_subsample as subsample {
      input:
        alignment_msa_fasta = mafft.aligned_sequences,
        sample_metadata_tsv = derived_cols.derived_metadata,
        build_name = build_name,
        builds_yaml = builds_yaml
    }
  }
  call utils.fasta_to_ids {
    input:
      sequences_fasta = select_first([subsample.subsampled_msa, mafft.aligned_sequences])
  }
  call nextstrain.snp_sites {
    input:
      msa_fasta = select_first([subsample.subsampled_msa, mafft.aligned_sequences])
  }
  #### augur_from_msa
  call nextstrain.augur_mask_sites {
    input:
      sequences = select_first([subsample.subsampled_msa, mafft.aligned_sequences])
  }
  call nextstrain.draft_augur_tree {
    input:
     msa_or_vcf = augur_mask_sites.masked_sequences
  }
  call phylo.snp_dists {
    input:
      cluster_name = build_name,
      alignment = select_first([subsample.subsampled_msa, mafft.aligned_sequences])
  }
  call versioning.version_capture{
    input:
  }
  output {
    # Version Capture
    String TheiaCoV_Augur_DistanceTree_version = version_capture.phvg_version
    String TheiaCoV_Augur_DistanceTree_analysis_date = version_capture.date
    # Tree, Intermediates, and Metadata
    File combined_assemblies = filter_sequences_by_length.filtered_fasta
    File multiple_alignment = mafft.aligned_sequences
    File unmasked_snps = snp_sites.snps_vcf
    File masked_alignment = augur_mask_sites.masked_sequences
    File metadata_merged = derived_cols.derived_metadata
    File keep_list = fasta_to_ids.ids_txt
    File mafft_alignment = select_first([subsample.subsampled_msa, mafft.aligned_sequences])
    File distance_tree = draft_augur_tree.aligned_tree
    # SNP Matrix
    File snp_matrix = snp_dists.snp_matrix
  }
}