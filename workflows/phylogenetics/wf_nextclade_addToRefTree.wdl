version 1.0

import "../../tasks/taxon_id/task_nextclade_add_ref.wdl" as nextclade_analysis

workflow nextclade_addToRefTree {
    meta {
      description: "Nextclade workflow that adds samples to a curated JSON tree from Augur."
    }
    input {
      File assembly_fastas
      File? root_sequence_fasta
      File? gene_annotations_gff
      File? reference_tree_json
      File? qc_config_json
      File? pcr_primers_csv
      File? virus_properties
      String dataset_name
      String? dataset_reference
      String? dataset_tag
    }
    call nextclade_analysis.nextclade { # nextclade analysis
      input:
        genome_fasta = assembly_fastas,
        root_sequence = root_sequence_fasta,
        reference_tree_json = reference_tree_json,
        qc_config_json = qc_config_json,
        gene_annotations_gff = gene_annotations_gff,
        pcr_primers_csv = pcr_primers_csv,
        virus_properties = virus_properties,
        dataset_name = dataset_name,
        dataset_reference = dataset_reference,
        dataset_tag = dataset_tag
    }
    output {
      String treeUpdate_nextclade_version = select_first([nextclade.nextclade_version, ""])
      File treeUpdate_nextclade_json = select_first([nextclade.nextclade_json, ""])
      File treeUpdate_auspice_json = select_first([nextclade.auspice_json, ""])
      File treeUpdate_nextclade_tsv = select_first([nextclade.nextclade_tsv, ""])
      String treeUpdate_nextclade_docker = select_first([nextclade.nextclade_docker, ""])
    }
}