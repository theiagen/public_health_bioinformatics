version 1.0

task set_organism_defaults_flu {
  input {
    String flu_segment
    String flu_subtype
    # nextclade references
    String ref_flu_h1n1_ha = "gs://theiagen-public-files-rp/terra/flu-references/reference_h1n1pdm_ha.fasta"
    String ref_flu_h1n1_na = "gs://theiagen-public-files-rp/terra/flu-references/reference_h1n1pdm_na.fasta"
    String ref_flu_h3n2_ha = "gs://theiagen-public-files-rp/terra/flu-references/reference_h3n2_ha.fasta"
    String ref_flu_h3n2_na = "gs://theiagen-public-files-rp/terra/flu-references/reference_h3n2_na.fasta"
    String ref_flu_yam_ha = "gs://theiagen-public-files-rp/terra/flu-references/reference_yam_ha.fasta"
    String ref_flu_yam_na = "gs://theiagen-public-files-rp/terra/flu-references/reference_yam_na.fasta"
    String ref_flu_vic_ha = "gs://theiagen-public-files-rp/terra/flu-references/reference_vic_ha.fasta"
    String ref_flu_vic_na = "gs://theiagen-public-files-rp/terra/flu-references/reference_vic_na.fasta"

    String nextclade_flu_h1n1_ha_tag = "2023-08-10T12:00:00Z"
    String nextclade_flu_h1n1_na_tag = "2023-08-10T12:00:00Z"
    String nextclade_flu_h3n2_ha_tag = "2023-08-10T12:00:00Z"
    String nextclade_flu_h3n2_na_tag = "2023-08-10T12:00:00Z"
    String nextclade_flu_vic_ha_tag = "2023-08-10T12:00:00Z"
    String nextclade_flu_vic_na_tag = "2023-08-10T12:00:00Z"

    Int genome_len = 13000
    Boolean trim_primer = false
  }
  command <<<
    echo "gs://theiagen-public-files/terra/theiacov-files/empty.fasta" > REFERENCE
    echo "NA" > NEXTCLADE_DATASET_TAG
    echo "NA" > NEXTCLADE_DATASET_NAME
    echo "NA" > NEXTCLADE_REFERENCE

    if [ "~{flu_segment}" == "HA" ]; then
      if [ "~{flu_subtype}" == "H1N1" ]; then
        echo "~{ref_flu_h1n1_ha}" > REFERENCE
        echo "~{nextclade_flu_h1n1_ha_tag}" > NEXTCLADE_DATASET_TAG
        echo "flu_h1n1pdm_ha" > NEXTCLADE_DATASET_NAME
        echo "MW626062" > NEXTCLADE_REFERENCE
      elif [ "~{flu_subtype}" == "H3N2" ]; then
        echo "~{ref_flu_h3n2_ha}" > REFERENCE
        echo "~{nextclade_flu_h3n2_ha_tag}" > NEXTCLADE_DATASET_TAG
        echo "flu_h3n2_ha" > NEXTCLADE_DATASET_NAME
        echo "EPI1857216" > NEXTCLADE_REFERENCE
      elif [ "~{flu_subtype}" == "Victoria" ]; then
        echo "~{ref_flu_vic_ha}" > REFERENCE
        echo "~{nextclade_flu_vic_ha_tag}" > NEXTCLADE_DATASET_TAG
        echo "flu_vic_ha" > NEXTCLADE_DATASET_NAME
        echo "KX058884" > NEXTCLADE_REFERENCE
      elif [ "~{flu_subtype}" == "Yamagata" ]; then
        echo "~{ref_flu_yam_ha}" > REFERENCE
        echo "flu_yam_ha" > NEXTCLADE_DATASET_NAME
        echo "JN993010" > NEXTCLADE_REFERENCE
      fi
    elif [ "~{flu_segment}" == "NA" ]; then
      if [ "~{flu_subtype}" == "H1N1" ]; then
        echo "~{ref_flu_h1n1_na}" > REFERENCE
        echo "~{nextclade_flu_h1n1_na_tag}" > NEXTCLADE_DATASET_TAG
        echo "flu_h1n1pdm_na" > NEXTCLADE_DATASET_NAME
        echo "MW626056" > NEXTCLADE_REFERENCE
      elif [ "~{flu_subtype}" == "H3N2" ]; then
        echo "~{ref_flu_h3n2_na}" > REFERENCE
        echo "~{nextclade_flu_h3n2_na_tag}" > NEXTCLADE_DATASET_TAG
        echo "flu_h3n2_na" > NEXTCLADE_DATASET_NAME
        echo "EPI1857215" > NEXTCLADE_REFERENCE
      elif [ "~{flu_subtype}" == "Victoria" ]; then
        echo "~{ref_flu_vic_na}" > REFERENCE
        echo "~{nextclade_flu_vic_na_tag}" > NEXTCLADE_DATASET_TAG
        echo "flu_vic_na" > NEXTCLADE_DATASET_NAME
        echo "CY073894" > NEXTCLADE_REFERENCE
      elif [ "~{flu_subtype}" == "Yamagata" ]; then
        echo "~{ref_flu_yam_na}" > REFERENCE
      fi
    fi
  >>>
  output {
    String reference = read_string("REFERENCE")
    String nextclade_dataset_tag = read_string("NEXTCLADE_DATASET_TAG")
    String nextclade_dataset_name = read_string("NEXTCLADE_DATASET_NAME")
    String nextclade_reference = read_string("NEXTCLADE_REFERENCE")
    Int genome_length = genome_len
    Boolean trim_primers = trim_primer
  }
  runtime {
    memory: "1 GB"
    cpu: 1
    docker: "ubuntu:jammy"
    disks: "local-disk 10 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2" 
  }
}