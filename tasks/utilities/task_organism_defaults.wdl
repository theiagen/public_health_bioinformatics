version 1.0

task set_organism_defaults_sc2 {
  input {
    String reference_genome = "gs://theiagen-public-files-rp/terra/augur-sars-cov-2-references/MN908947.fasta"
    String nextclade_ds_tag = "2023-08-17T12:00:00Z"
    String nextclade_ref = "MN908947"
    String nextclade_ds_name = "sars-cov-2"
    Int genome_len = 29903
    Int vadr_max_length = 30000
    String vadr_options = "--noseqnamemax --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn --out_allfasta"
  }
  command <<<
    echo "setting defaults for sc2"
  >>>
  output {
    String reference = reference_genome
    String nextclade_dataset_tag = nextclade_ds_tag
    String nextclade_reference = nextclade_ref
    String nextclade_dataset_name = nextclade_ds_name
    Int genome_length = genome_len
    Int vadr_maxlen = vadr_max_length
    String vadr_opts = vadr_options
  }
  runtime {
    memory: "1 GB"
    cpu: 1
    docker: "ubuntu:jammy"
    disks: "local-disk 10 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2" 
  }
}

task set_organism_defaults_mpox {
  input {
    String reference_genome = "gs://theiagen-public-files/terra/mpxv-files/MPXV.MT903345.reference.fasta"
    String nextclade_ds_tag = "2023-08-01T12:00:00Z"
    String nextclade_ref = "pseudo_ON563414"
    String nextclade_ds_name = "hMPXV_B1"
    String target_org = "Monkeypox virus"
    String primer_bed_file = "gs://theiagen-public-files/terra/mpxv-files/MPXV.primer.bed"
    String reference_gff_file = "gs://theiagen-public-files/terra/mpxv-files/Mpox-MT903345.1.reference.gff3"
    String vadr_options = "--glsearch -s -r --nomisc --mkey mpxv --r_lowsimok --r_lowsimxd 100 --r_lowsimxl 2000 --alt_pass discontn,dupregin --out_allfasta --minimap2 --s_overhang 150"
    Int vadr_max_length = 210000
    Int genome_len = 197200
  }
  command <<<
    echo "setting defaults for mpox"
  >>>
  output {
    String reference = reference_genome
    String nextclade_dataset_tag = nextclade_ds_tag
    String nextclade_reference = nextclade_ref
    String nextclade_dataset_name = nextclade_ds_name
    String target_organism = target_org
    String primer_bed = primer_bed_file
    String reference_gff = reference_gff_file
    String vadr_opts = vadr_options
    Int vadr_maxlen = vadr_max_length
    Int genome_length = genome_len
  }
  runtime {
    memory: "1 GB"
    cpu: 1
    docker: "ubuntu:jammy"
    disks: "local-disk 10 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2" 
  }
}
 
task set_organism_defaults_wnv {
  input {
    String reference_genome = "gs://theiagen-public-files/terra/theiacov-files/WNV/NC_009942.1_wnv_L1.fasta"
    String target_org = "West Nile virus"
    String primer_bed_file = "gs://theiagen-public-files/terra/theiacov-files/WNV/WNV-L1_primer.bed"
    Int genome_len = 11000
    String vadr_options = "--mkey flavi --mdir /opt/vadr/vadr-models-flavi/ --nomisc --noprotid --out_allfasta"    
    Int vadr_max_length = 11000
  }
  command <<<
    echo "setting defaults for wnv"
  >>>
  output {
    String reference = reference_genome
    String target_organism = target_org
    String primer_bed = primer_bed_file    
    Int genome_length = genome_len
    String vadr_opts = vadr_options
    Int vadr_maxlen = vadr_max_length
  }
  runtime {
    memory: "1 GB"
    cpu: 1
    docker: "ubuntu:jammy"
    disks: "local-disk 10 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2" 
  }
}

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
   
task set_organism_defaults_rsv_a {
  input {
    String reference_genome = "gs://theiagen-public-files-rp/terra/rsv_references/reference_rsv_a.fasta"
    String nextclade_ds_tag = "2023-02-03T12:00:00Z"
    String nextclade_ref = "EPI_ISL_412866"
    String nextclade_ds_name = "rsv_a"
    Int genome_len = 16000
    String vadr_options = "-r --mkey rsv --xnocomp"
    Int vadr_max_length = 15500
  }
  command <<<
    echo "setting defaults for rsva"
  >>>
  output {
    String reference = reference_genome
    String nextclade_dataset_tag = nextclade_ds_tag
    String nextclade_reference = nextclade_ref
    String nextclade_dataset_name = nextclade_ds_name
    Int genome_length = genome_len
    String vadr_opts = vadr_options
    Int vadr_maxlen = vadr_max_length
  }
  runtime {
    memory: "1 GB"
    cpu: 1
    docker: "ubuntu:jammy"
    disks: "local-disk 10 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2" 
  }
}

task set_organism_defaults_rsv_b {
  input {
    String reference_genome = "gs://theiagen-public-files-rp/terra/rsv_references/reference_rsv_b.fasta"
    String nextclade_ds_tag = "2023-02-03T12:00:00Z"
    String nextclade_ref = "EPI_ISL_1653999"
    String nextclade_ds_name = "rsv_b"
    Int genome_len = 16000   
    String vadr_options = "-r --mkey rsv --xnocomp"
    Int vadr_max_length = 15500
  }
  command <<<
    echo "setting defaults for rsvb"
  >>>
  output {
    String reference = reference_genome
    String nextclade_dataset_tag = nextclade_ds_tag
    String nextclade_reference = nextclade_ref
    String nextclade_dataset_name = nextclade_ds_name
    Int genome_length = genome_len
    String vadr_opts = vadr_options
    Int vadr_maxlen = vadr_max_length
  }
  runtime {
    memory: "1 GB"
    cpu: 1
    docker: "ubuntu:jammy"
    disks: "local-disk 10 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2" 
  }
}