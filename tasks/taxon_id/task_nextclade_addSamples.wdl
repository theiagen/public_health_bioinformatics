version 1.0

task nextclade {
    meta {
      description: "Nextclade classification of one sample. Leaving optional inputs unspecified will use SARS-CoV-2 defaults."
    }
    input {
      File genome_fasta
      File? root_sequence
      File? reference_tree_json
      File? qc_config_json
      File? gene_annotations_gff
      File? pcr_primers_csv
      File? virus_properties
      String docker = "us-docker.pkg.dev/general-theiagen/nextstrain/nextclade:2.14.0"
      String dataset_name
      String? dataset_reference
      String? dataset_tag
      Int disk_size = 50
    }
    String basename = basename(genome_fasta, ".fasta")
    command <<<
        NEXTCLADE_VERSION="$(nextclade --version)"
        echo $NEXTCLADE_VERSION > NEXTCLADE_VERSION

        nextclade dataset get \
          --name="~{dataset_name}" \
          ~{"--reference " + dataset_reference} \
          ~{"--tag " + dataset_tag} \
          -o nextclade_dataset_dir \
          --verbose
        
        # If no referece sequence is provided, use the reference tree from the dataset
        if [[ ! -z "~{reference_tree_json}" ]]; then
          reference_tree_json=nextclade_dataset_dir/tree.json
        else
          reference_tree_json="~{reference_tree_json}"
        fi

        set -e
        nextclade run \
            --input-dataset=nextclade_dataset_dir/ \
            ~{"--input-root-seq " + root_sequence} \
            --input-tree ${reference_tree_json} \
            ~{"--input-qc-config " + qc_config_json} \
            ~{"--input-gene-map " + gene_annotations_gff} \
            ~{"--input-pcr-primers " + pcr_primers_csv} \
            ~{"--input-virus-properties " + virus_properties}  \
            --output-json "~{basename}".nextclade.json \
            --output-tsv  "~{basename}".nextclade.tsv \
            --output-tree "~{basename}".nextclade.auspice.json \
            --output-all=. \
            "~{genome_fasta}"
    >>>
    runtime {
      docker: "~{docker}"
      memory: "8 GB"
      cpu: 2
      disks:  "local-disk " + disk_size + " SSD"
      disk: disk_size + " GB" # TES
      dx_instance_type: "mem1_ssd1_v2_x2"
      maxRetries: 3 
    }
    output {
      String nextclade_version = read_string("NEXTCLADE_VERSION")
      File nextclade_json = "~{basename}.nextclade.json"
      File auspice_json = "~{basename}.nextclade.auspice.json"
      File nextclade_tsv = "~{basename}.nextclade.tsv"
      String nextclade_docker = docker
    }
}