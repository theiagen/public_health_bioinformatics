version 1.0

import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade_analysis

workflow nextclade_addToRefTree {
    meta {
      description: "Nextclade workflow that adds samples to a curated JSON tree from Augur."
    }
    input {
      File assembly_fastas
      String build_name
      File root_sequence_fasta
      #File? gene_annotations_gff
      File reference_tree_json
      File? qc_config_json
      File? pcr_primers_csv
      File? virus_properties
      String docker = "nextstrain/nextclade:2.13.0"
      String dataset_name = "MPXV"
      String dataset_reference = "ancestral"
      String dataset_tag = "2023-01-26T12:00:00Z"
    }
    call nextclade_analysis.nextclade { # nextclade analysis
      input:
        genome_fasta = assembly_fastas,
        root_sequence = root_sequence_fasta,
        auspice_reference_tree_json = reference_tree_json,
        qc_config_json = qc_config_json,
        #gene_annotations_json = gene_annotations_gff,
        pcr_primers_csv = pcr_primers_csv,
        virus_properties = virus_properties,
        docker = docker,
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