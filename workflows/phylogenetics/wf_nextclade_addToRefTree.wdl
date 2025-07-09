version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade_analysis

workflow nextclade_addToRefTree {
    meta {
      description: "Nextclade workflow that adds samples to a curated JSON tree from Augur."
    }
    input {
      Array[File] assembly_fastas
      File? input_ref
      File? gene_annotations_gff
      File? reference_tree_json
      File? nextclade_pathogen_json
      String nextclade_dataset_name
      String? dataset_tag
    }
    call nextclade_analysis.nextclade_add_ref { # nextclade analysis
      input:
        genome_fastas = assembly_fastas,
        reference_tree_json = reference_tree_json,
        nextclade_pathogen_json = nextclade_pathogen_json,
        gene_annotations_gff = gene_annotations_gff,
        input_ref = input_ref,
        dataset_name = nextclade_dataset_name,
        dataset_tag = dataset_tag
    }
    call versioning.version_capture {
      input:
    }
    output {
      String treeUpdate_nextclade_version = select_first([nextclade_add_ref.nextclade_version, ""])
      File treeUpdate_nextclade_json = select_first([nextclade_add_ref.nextclade_json, ""])
      File treeUpdate_auspice_json = select_first([nextclade_add_ref.auspice_json, ""])
      File treeUpdate_nextclade_tsv = select_first([nextclade_add_ref.nextclade_tsv, ""])
      String treeUpdate_nextclade_docker = select_first([nextclade_add_ref.nextclade_docker, ""])
      # Version Capture
      String samples_to_ref_tree_version = version_capture.phb_version
      String samples_to_ref_tree_analysis_date = version_capture.date
    }
}
