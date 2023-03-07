version 1.0

import "wf_sarscov2_nextstrain_modified.wdl" as augur
import "../../tasks/phylogenetic_inference/task_snp_dists.wdl" as snp_dists_task
import "../../tasks/task_versioning.wdl" as versioning

workflow theiacov_augur_run {
  meta {
    description: "Workflow for SC2 cluster investigations. TheiaCoV_Augur_Run will run Augur without a subsampling module using a modified version of The Broad Institute's sarscov2_nextstrain WDL workflow to create an Auspice JSON file; output from the modified sarscov2_nextstrain workflow will also be used to infer SNP distances"
    author: "Kevin G Libuit"
    email:  "kevin.libuit@theiagen.com"
  }
  input {
    Array[File]+ assembly_fastas
    Array[File]+ sample_metadata_tsvs
    String build_name
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
  }
  call augur.sarscov2_nextstrain {
    input:
      assembly_fastas = assembly_fastas,
      sample_metadata_tsvs = sample_metadata_tsvs,
      build_name = build_name
  } 
  call snp_dists_task.snp_dists {
    input:
      cluster_name = build_name,
      alignment = sarscov2_nextstrain.mafft_alignment
  }
  call versioning.version_capture {
    input:
  }
  output {
    # Version Capture
    String theiacov_augur_run_version = version_capture.phb_version
    String theiacov_augur_run_analysis_date = version_capture.date
    # Augur outputs
    File combined_assemblies = sarscov2_nextstrain.combined_assemblies
    File MAFFT_alignment = sarscov2_nextstrain.mafft_alignment
    File unmasked_snps = sarscov2_nextstrain.unmasked_snps
    File metadata_merged = sarscov2_nextstrain.metadata_merged
    File keep_list = sarscov2_nextstrain.keep_list
    File distance_tree = sarscov2_nextstrain.ml_tree
    File time_tree = sarscov2_nextstrain.time_tree
    File auspice_input_json = sarscov2_nextstrain.auspice_input_json
    # SNP Matrix
    File snp_matrix = snp_dists.snp_matrix
  }
}