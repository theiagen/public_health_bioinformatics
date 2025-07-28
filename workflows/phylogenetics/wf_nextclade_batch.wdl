version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade_analysis

workflow nextclade_batch {
    meta {
      description: "Nextclade workflow that adds samples to a curated JSON tree from Augur."
    }
    input {
      Array[File] assembly_fastas
      File? input_ref
      File? gene_annotations_gff
      File? reference_tree_json
      File? pathogen_json
      String dataset_name
      String? dataset_tag
    }
    call nextclade_analysis.nextclade_v3_set { # nextclade analysis
      input:
        genome_fastas = assembly_fastas,
        reference_tree_json = reference_tree_json,
        pathogen_json = pathogen_json,
        gene_annotations_gff = gene_annotations_gff,
        input_ref = input_ref,
        dataset_name = dataset_name,
        dataset_tag = dataset_tag
    }
    call versioning.version_capture {
      input:
    }
    output {
      String nextclade_batch_nextclade_version = select_first([nextclade_v3_set.nextclade_version, ""])
      File nextclade_batch_nextclade_json = select_first([nextclade_v3_set.nextclade_json, ""])
      File nextclade_batch_auspice_json = select_first([nextclade_v3_set.auspice_json, ""])
      File nextclade_batch_nextclade_tsv = select_first([nextclade_v3_set.nextclade_tsv, ""])
      String nextclade_batch_nextclade_docker = select_first([nextclade_v3_set.nextclade_docker, ""])
      # Version Capture
      String nextclade_batch_version = version_capture.phb_version
      String nextclade_batch_analysis_date = version_capture.date
    }
}
