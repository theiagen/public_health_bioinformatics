version 1.0

workflow organism_parameters {
  meta {
    description: "Organizes all organism defaults into a single location for easier maintanence"
  }
  input {
    String organism
    
    # hiv information
    String hiv_primer_version = "v1"

    # flu information
    String flu_segment = "N/A"
    String flu_subtype = "N/A"

    # sequencing & reference information
    File? primer_bed_file
    File? reference_gff_file
    File? reference_genome
    File? gene_locations_bed_file
    Int? genome_length_input

    # set default nextclade information as NA
    String? nextclade_dataset_tag_input
    String? nextclade_dataset_name_input

    # vadr parameters
    Int? vadr_max_length
    String? vadr_options
    Int? vadr_mem

    # pangolin parameters
    String? pangolin_docker_image

    # kraken parameters
    String? kraken_target_organism_input
  }
  if (organism == "sars-cov-2" || organism == "SARS-CoV-2") {
    String sc2_org_name = "sars-cov-2"
    String sc2_reference_genome = "gs://theiagen-public-files-rp/terra/augur-sars-cov-2-references/MN908947.fasta"
    String sc2_gene_locations_bed = "gs://theiagen-public-files-rp/terra/sars-cov-2-files/sc2_gene_locations.bed"
    String sc2_nextclade_ds_tag = "2024-04-15--15-08-22Z"
    String sc2_nextclade_ds_name = "nextstrain/sars-cov-2/wuhan-hu-1/orfs"
    String sc2_pangolin_docker = "us-docker.pkg.dev/general-theiagen/staphb/pangolin:4.3.1-pdata-1.26"
    Int sc2_genome_len = 29903
    Int sc2_vadr_max_length = 30000
    String sc2_vadr_options = "--noseqnamemax --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn --out_allfasta"
    Int sc2_vadr_memory = 8
  }
  if (organism == "MPXV" || organism == "mpox" || organism == "monkeypox" || organism == "Monkeypox virus" || organism == "Mpox") {
    String mpox_org_name = "MPXV"
    String mpox_reference_genome = "gs://theiagen-public-files/terra/mpxv-files/MPXV.MT903345.reference.fasta"
    String mpox_gene_locations_bed = "gs://theiagen-public-files/terra/mpxv-files/mpox_gene_locations.bed"
    String mpox_nextclade_ds_tag = "2024-01-16--20-31-02Z"
    String mpox_nextclade_ds_name = "nextstrain/mpox/lineage-b.1"
    String mpox_kraken_target_organism = "Monkeypox virus"
    String mpox_primer_bed_file = "gs://theiagen-public-files/terra/mpxv-files/MPXV.primer.bed"
    String mpox_reference_gff_file = "gs://theiagen-public-files/terra/mpxv-files/Mpox-MT903345.1.reference.gff3"
    String mpox_vadr_options = "--glsearch -s -r --nomisc --mkey mpxv --r_lowsimok --r_lowsimxd 100 --r_lowsimxl 2000 --alt_pass discontn,dupregin --out_allfasta --minimap2 --s_overhang 150"
    Int mpox_vadr_max_length = 210000
    Int mpox_vadr_memory = 8
    Int mpox_genome_len = 197200
  }  
  if (organism == "WNV" || organism == "wnv" || organism == "West Nile virus") {
    String wnv_org_name = "WNV"
    String wnv_reference_genome = "gs://theiagen-public-files/terra/theiacov-files/WNV/NC_009942.1_wnv_L1.fasta"
    String wnv_kraken_target_organism = "West Nile virus"
    String wnv_primer_bed_file = "gs://theiagen-public-files/terra/theiacov-files/WNV/WNV-L1_primer.bed"
    Int wnv_genome_len = 11000
    String wnv_vadr_options = "--mkey flavi --mdir /opt/vadr/vadr-models-flavi/ --nomisc --noprotid --out_allfasta"    
    Int wnv_vadr_max_length = 11000
    Int wnv_vadr_memory = 8
    String wnv_nextclade_ds_tag = "NA"
    String wnv_nextclade_ds_name = "NA"
  }
  if (organism == "flu" || organism == "influenza" || organism == "Flu" || organism == "Influenza") {
    String flu_org_name = "flu"
    Int flu_genome_len = 13500

    # vadr options are dummy options for flu right now
    String flu_vadr_options = "--atgonly --xnocomp --nomisc --alt_fail extrant5,extrant3 --mkey flu"
    Int flu_vadr_max_length = 13500
    Int flu_vadr_memory = 8

    # setting nextclade parameters
    if (flu_segment == "HA") {
      if (flu_subtype == "H1N1") {
        String h1n1_ha_reference = "gs://theiagen-public-files-rp/terra/flu-references/reference_h1n1pdm_ha.fasta"
        String h1n1_ha_nextclade_ds_tag = "2024-01-16--20-31-02Z"
        String h1n1_ha_nextclade_ds_name = "nextstrain/flu/h1n1pdm/ha/MW626062"
      }
      if (flu_subtype == "H3N2") {
        String h3n2_ha_reference = "gs://theiagen-public-files-rp/terra/flu-references/reference_h3n2_ha.fasta"
        String h3n2_ha_nextclade_ds_tag = "2024-02-22--16-12-03Z"
        String h3n2_ha_nextclade_ds_name = "nextstrain/flu/h3n2/ha/EPI1857216"
      }
      if (flu_subtype == "Victoria") {
        String vic_ha_reference = "gs://theiagen-public-files-rp/terra/flu-references/reference_vic_ha.fasta"
        String vic_ha_nextclade_ds_tag = "2024-01-16--20-31-02Z"
        String vic_ha_nextclade_ds_name = "nextstrain/flu/vic/ha/KX058884"
      }
      if (flu_subtype == "Yamagata") {
        String yam_ha_reference = "gs://theiagen-public-files-rp/terra/flu-references/reference_yam_ha.fasta"
        String yam_ha_nextclade_ds_tag = "2024-01-30--16-34-55Z"
        String yam_ha_nextclade_ds_name = "nextstrain/flu/yam/ha/JN993010"
      }
    }
    if (flu_segment == "NA") {
      if (flu_subtype == "H1N1") {
        String h1n1_na_reference = "gs://theiagen-public-files-rp/terra/flu-references/reference_h1n1pdm_na.fasta"
        String h1n1_na_nextclade_ds_tag = "2024-01-16--20-31-02Z"
        String h1n1_na_nextclade_ds_name = "nextstrain/flu/h1n1pdm/na/MW626056"
      }
      if (flu_subtype == "H3N2") {
        String h3n2_na_reference = "gs://theiagen-public-files-rp/terra/flu-references/reference_h3n2_na.fasta"
        String h3n2_na_nextclade_ds_tag = "2024-01-16--20-31-02Z"
        String h3n2_na_nextclade_ds_name = "nextstrain/flu/h3n2/na/EPI1857215"
      }
      if (flu_subtype == "Victoria") {
        String vic_na_reference = "gs://theiagen-public-files-rp/terra/flu-references/reference_vic_na.fasta"
        String vic_na_nextclade_ds_tag = "2024-01-16--20-31-02Z"
        String vic_na_nextclade_ds_name = "nextstrain/flu/vic/na/CY073894"
      }
      if (flu_subtype == "Yamagata") {
        String yam_na_reference = "gs://theiagen-public-files-rp/terra/flu-references/reference_yam_na.fasta"
        String yam_na_nextclade_ds_tag = "NA"
        String yam_na_nextclade_ds_name = "NA"
      }
    }
  }
  if (organism == "rsv_a" || organism == "rsv-a" || organism == "RSV-A" || organism == "RSV_A") {
    String rsv_a_org_name = "rsv_a"
    String rsv_a_reference_genome = "gs://theiagen-public-files-rp/terra/rsv_references/reference_rsv_a.fasta"
    String rsv_a_nextclade_ds_tag = "2024-01-29--10-29-43Z"
    String rsv_a_nextclade_ds_name = "nextstrain/rsv/a/EPI_ISL_412866"
    Int rsv_a_genome_len = 16000
    String rsv_a_vadr_options = "-r --mkey rsv --xnocomp"
    Int rsv_a_vadr_max_length = 15500
    Int rsv_a_vadr_memory = 32
  }
  if (organism == "rsv_b" || organism == "rsv-b" || organism == "RSV-B" || organism == "RSV_B") {
    String rsv_b_org_name = "rsv_b"
    String rsv_b_reference_genome = "gs://theiagen-public-files-rp/terra/rsv_references/reference_rsv_b.fasta"
    String rsv_b_nextclade_ds_tag = "2024-01-29--10-29-43Z"
    String rsv_b_nextclade_ds_name = "nextstrain/rsv/b/EPI_ISL_1653999"
    Int rsv_b_genome_len = 16000   
    String rsv_b_vadr_options = "-r --mkey rsv --xnocomp"
    Int rsv_b_vadr_max_length = 15500
    Int rsv_b_vadr_memory = 32
  }
  if (organism == "HIV" && hiv_primer_version == "v1") {
    String hiv_v1_org_name = "HIV"
    String hiv_v1_reference_genome = "gs://theiagen-public-files/terra/hivgc-files/NC_001802.1.fasta"
    String hiv_v1_reference_gff = "gs://theiagen-public-files/terra/hivgc-files/NC_001802.1.gff3"
    String hiv_v1_primer_bed = "gs://theiagen-public-files/terra/hivgc-files/HIV-1_v1.0.primer.hyphen.bed"
    String hiv_v1_target_organism = "Human immunodeficiency virus 1"
    Int hiv_v1_genome_len = 9181 
  }
  if (organism == "HIV" && hiv_primer_version == "v2") {
    String hiv_v2_org_name = "HIV"
    String hiv_v2_reference_genome = "gs://theiagen-public-files/terra/hivgc-files/AY228557.1.headerchanged.fasta"
    String hiv_v2_reference_gff = "gs://theiagen-public-files/terra/hivgc-files/AY228557.1.gff3"
    String hiv_v2_primer_bed = "gs://theiagen-public-files/terra/hivgc-files/HIV-1_v2.0.primer.hyphen400.1.bed"
    String hiv_v2_target_organism = "Human immunodeficiency virus 1"
    Int hiv_v2_genome_len = 9840
  }
  output {
    # standardized organism flag
    String standardized_organism = select_first([sc2_org_name, mpox_org_name, wnv_org_name, flu_org_name, rsv_a_org_name, rsv_b_org_name, hiv_v1_org_name, hiv_v2_org_name])
    # reference genome and sequencing information
    File reference = select_first([reference_genome, sc2_reference_genome, mpox_reference_genome, wnv_reference_genome, h1n1_ha_reference, h3n2_ha_reference, vic_ha_reference, yam_ha_reference, h1n1_na_reference, h3n2_na_reference, vic_na_reference, yam_na_reference, 
    rsv_a_reference_genome, rsv_b_reference_genome, hiv_v1_reference_genome, hiv_v2_reference_genome, "gs://theiagen-public-files/terra/theiacov-files/empty.fasta"])
    File gene_locations_bed = select_first([gene_locations_bed_file, sc2_gene_locations_bed, mpox_gene_locations_bed, "gs://theiagen-public-files/terra/theiacov-files/empty.bed"])
    File primer_bed = select_first([primer_bed_file, mpox_primer_bed_file, wnv_primer_bed_file, hiv_v1_primer_bed, hiv_v2_primer_bed, "gs://theiagen-public-files/terra/theiacov-files/empty.bed"])
    File reference_gff = select_first([reference_gff_file, mpox_reference_gff_file, hiv_v1_reference_gff, hiv_v2_reference_gff, "gs://theiagen-public-files/terra/theiacov-files/empty.gff3"])
    Int genome_length = select_first([genome_length_input, sc2_genome_len, mpox_genome_len, wnv_genome_len, flu_genome_len, rsv_a_genome_len, rsv_b_genome_len, hiv_v1_genome_len, hiv_v2_genome_len])
    # nextclade information
    String nextclade_dataset_tag = select_first([nextclade_dataset_tag_input, sc2_nextclade_ds_tag, mpox_nextclade_ds_tag, wnv_nextclade_ds_tag, h1n1_ha_nextclade_ds_tag, h3n2_ha_nextclade_ds_tag, vic_ha_nextclade_ds_tag, yam_ha_nextclade_ds_tag, h1n1_na_nextclade_ds_tag, h3n2_na_nextclade_ds_tag, vic_na_nextclade_ds_tag, yam_na_nextclade_ds_tag, rsv_a_nextclade_ds_tag, rsv_b_nextclade_ds_tag, "NA"])
    String nextclade_dataset_name = select_first([nextclade_dataset_name_input, sc2_nextclade_ds_name, mpox_nextclade_ds_name, wnv_nextclade_ds_name, h1n1_ha_nextclade_ds_name, h3n2_ha_nextclade_ds_name, vic_ha_nextclade_ds_name, yam_ha_nextclade_ds_name, h1n1_na_nextclade_ds_name, h3n2_na_nextclade_ds_name, vic_na_nextclade_ds_name, yam_na_nextclade_ds_name, rsv_a_nextclade_ds_name, rsv_b_nextclade_ds_name, "NA"])
    # pangolin options
    String pangolin_docker = select_first([pangolin_docker_image, sc2_pangolin_docker, ""])
    # vadr options
    String vadr_opts = select_first([vadr_options, sc2_vadr_options, mpox_vadr_options, wnv_vadr_options, flu_vadr_options, rsv_a_vadr_options, rsv_b_vadr_options, "NA"])
    Int vadr_maxlength = select_first([vadr_max_length, sc2_vadr_max_length, mpox_vadr_max_length, wnv_vadr_max_length, flu_vadr_max_length, rsv_a_vadr_max_length, rsv_b_vadr_max_length, 0])
    Int vadr_memory = select_first([vadr_mem, sc2_vadr_memory, mpox_vadr_memory, wnv_vadr_memory, flu_vadr_memory, rsv_a_vadr_memory, rsv_b_vadr_memory, 0])
    # kraken options
    String kraken_target_organism = select_first([kraken_target_organism_input, mpox_kraken_target_organism, wnv_kraken_target_organism, hiv_v1_target_organism, hiv_v2_target_organism, ""])
  }
}