version 1.0

workflow augur_parameters {
  meta {
    description: "Organizes all augur organism defaults into a single location for easier maintanence"
  }
  input {
    String organism

    # flu information
    String flu_segment = "N/A"
    String flu_subtype = "N/A"
    String flu_genoflu_genotype = "N/A"

    File? reference_genome # this is for the most (all?) part fasta, so we may opt to rename to reference_fasta
    File? reference_genbank # for augur

    # augur parameters
    Int? min_num_unambig
    File? clades_tsv
    File? auspice_config
    Int? pivot_interval
    Float? min_date
    Float? narrow_bandwidth
    Float? proportion_wide
  }
  if (organism == "sars-cov-2" || organism == "SARS-CoV-2" || organism == "2697049" || organism == "3418604") {
    String sc2_org_name = "sars-cov-2"
    String sc2_reference_genome = "gs://theiagen-public-resources-rp/reference_data/viral/sars-cov-2/MN908947.fasta"

    File sc2_clades_tsv = "gs://theiagen-public-resources-rp/reference_data/viral/sars-cov-2/sc2_clades_20251008.tsv"
    File sc2_auspice_config = "gs://theiagen-public-resources-rp/reference_data/viral/sars-cov-2/sc2_auspice_config_20251030.json"
    File sc2_reference_genbank = "gs://theiagen-public-resources-rp/reference_data/viral/sars-cov-2/sc2_reference_seq_20251008.gb"
    Float sc2_min_date = 2020.0
    Int sc2_pivot_interval = 1
    String sc2_pivot_interval_units = "weeks"
    Float sc2_narrow_bandwidth = 0.05
    Float sc2_proportion_wide = 0.0
    Int sc2_min_num_unambig = 27000
  }
  if (organism == "MPXV" || organism == "mpox" || organism == "monkeypox" || organism == "Monkeypox virus" || organism == "Mpox" || organism == "10244") {
    String mpox_org_name = "MPXV"
    String mpox_reference_genome = "gs://theiagen-public-resources-rp/reference_data/viral/mpox/MPXV.MT903345.reference.fasta"

    File mpox_clades_tsv = "gs://theiagen-public-resources-rp/reference_data/viral/mpox/mpox_clades.tsv"
    File mpox_reference_gbk = "gs://theiagen-public-resources-rp/reference_data/viral/mpox/NC_063383.1_reference.gb"
    File mpox_auspice_config = "gs://theiagen-public-resources-rp/reference_data/viral/mpox/mpox_auspice_config_mpxv.json"
    Int mpox_min_num_unambig = 150000
    # inherited from flu defaults
    Float mpox_min_date = 2020.0
    Int mpox_pivot_interval = 1
    Float mpox_narrow_bandwidth = 0.1666667
    Float mpox_proportion_wide = 0.0
  }
  if (organism == "flu" || organism == "influenza" || organism == "Flu" || organism == "Influenza" || organism == "11320" || organism == "11309" || organism == "11308" || organism == "11520") {
    String flu_org_name = "flu"
    Int flu_genome_len = 13500

    # augur options for flu
    Int flu_min_num_unambig = 900
    Float flu_min_date = 2020.0
    Int flu_pivot_interval = 1
    Float flu_narrow_bandwidth = 0.1666667
    Float flu_proportion_wide = 0.0

    # setting nextclade and augur parameters
    if (flu_segment == "HA") {
      if (flu_subtype == "H1N1") {
        String h1n1_ha_reference = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h1n1pdm_ha.fasta"
        String h1n1_ha_reference_gbk = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h1n1pdm_ha.gb"
        String h1n1_ha_clades_tsv = "gs://theiagen-public-resources-rp/reference_data/viral/flu/clades_h1n1pdm_ha.tsv"
        String h1n1_ha_auspice_config = "gs://theiagen-public-resources-rp/reference_data/viral/flu/auspice_config_h1n1pdm.json"
      }
      if (flu_subtype == "H3N2") {
        String h3n2_ha_reference = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h3n2_ha.fasta"
        String h3n2_ha_reference_gbk = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h3n2_ha.gb"
        String h3n2_ha_clades_tsv = "gs://theiagen-public-resources-rp/reference_data/viral/flu/clades_h3n2_ha.tsv"
        String h3n2_ha_auspice_config = "gs://theiagen-public-resources-rp/reference_data/viral/flu/auspice_config_h3n2.json"
      }
      if (flu_subtype == "Victoria") {
        String vic_ha_reference = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_vic_ha.fasta"
        String vic_ha_reference_gbk = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_vic_ha.gb"
        String vic_ha_clades_tsv = "gs://theiagen-public-resources-rp/reference_data/viral/flu/clades_vic_ha.tsv"
        String vic_ha_auspice_config = "gs://theiagen-public-resources-rp/reference_data/viral/flu/auspice_config_vic.json"
      }
      if (flu_subtype == "Yamagata") {
        String yam_ha_reference = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_yam_ha.fasta"
        String yam_ha_reference_gbk = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_yam_ha.gb"
        String yam_ha_clades_tsv = "gs://theiagen-public-resources-rp/reference_data/viral/flu/clades_yam_ha.tsv"
        String yam_ha_auspice_config = "gs://theiagen-public-resources-rp/reference_data/viral/flu/auspice_config_yam_20251030.json"
      }
      if (sub(flu_subtype, "^H5N.*$", "MATCH") == "MATCH" || flu_subtype == "H5") {
        # H5N1 is a special case where the dataset used is the h5nx all clades dataset
        String h5n1_ha_reference = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h5n1_ha.fasta"
        String h5n1_ha_reference_gbk = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h5n1_ha.gb"
        String h5n1_ha_auspice_config = "gs://theiagen-public-resources-rp/reference_data/viral/flu/auspice_config_h5n1_20251030.json"
      }
    }
    if (flu_segment == "NA") {
      if (flu_subtype == "H1N1") {
        String h1n1_na_reference = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h1n1pdm_na.fasta"
        String h1n1_na_reference_gbk = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h1n1pdm_na.gb"
        String h1n1_na_auspice_config = "gs://theiagen-public-resources-rp/reference_data/viral/flu/auspice_config_h1n1pdm_20251030.json"
      }
      if (flu_subtype == "H3N2") {
        String h3n2_na_reference = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h3n2_na.fasta"
        String h3n2_na_reference_gbk = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h3n2_na.gb"
        String h3n2_na_auspice_config = "gs://theiagen-public-resources-rp/reference_data/viral/flu/auspice_config_h3n2_20251030.json"
      }
      if (flu_subtype == "Victoria") {
        String vic_na_reference = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_vic_na.fasta"
        String vic_na_reference_gbk = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_yam_na.gb"
        String vic_na_auspice_config = "gs://theiagen-public-resources-rp/reference_data/viral/flu/auspice_config_vic_20251030.json"
      }
      if (flu_subtype == "Yamagata") {
        String yam_na_reference = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_yam_na.fasta"
        String yam_na_reference_gbk = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_yam_na.gb"
        String yam_na_auspice_config = "gs://theiagen-public-resources-rp/reference_data/viral/flu/auspice_config_yam_20251030.json"
      }
    }
  }
  if (organism == "rsv_a" || organism == "rsv-a" || organism == "RSV-A" || organism == "RSV_A" || organism == "208893") {
    String rsv_a_org_name = "rsv_a"
    String rsv_a_reference_genome = "gs://theiagen-public-resources-rp/reference_data/viral/rsv/reference_rsv_a.EPI_ISL_412866.fasta"

    # augur options for rsv-a
    File rsv_a_clades_tsv = "gs://theiagen-public-resources-rp/reference_data/viral/rsv/rsv_a_clades.tsv"
    File rsv_a_reference_gbk = "gs://theiagen-public-resources-rp/reference_data/viral/rsv/reference_rsv_a.gb"
    File rsv_a_auspice_config = "gs://theiagen-public-resources-rp/reference_data/viral/rsv/rsv_auspice_config_20251030.json"
    Int rsv_a_min_num_unambig = 10850 # using 70% of 15500
    # inherited from flu defaults
    Float rsv_a_min_date = 2020.0
    Int rsv_a_pivot_interval = 1
    Float rsv_a_narrow_bandwidth = 0.1666667
    Float rsv_a_proportion_wide = 0.0
  }
  if (organism == "rsv_b" || organism == "rsv-b" || organism == "RSV-B" || organism == "RSV_B" || organism == "208895") {
    String rsv_b_org_name = "rsv_b"
    String rsv_b_reference_genome = "gs://theiagen-public-resources-rp/reference_data/viral/rsv/reference_rsv_b.EPI_ISL_1653999.fasta"

    # augur options for rsv-b
    File rsv_b_clades_tsv = "gs://theiagen-public-resources-rp/reference_data/viral/rsv/rsv_b_clades.tsv"
    File rsv_b_reference_gbk = "gs://theiagen-public-resources-rp/reference_data/viral/rsv/reference_rsv_b.gb"
    File rsv_b_auspice_config = "gs://theiagen-public-resources-rp/reference_data/viral/rsv/rsv_auspice_config_20251030.json"
    Int rsv_b_min_num_unambig = 10850 # using 70% of 15500
    # inherited from flu defaults
    Float rsv_b_min_date = 2020.0
    Int rsv_b_pivot_interval = 1
    Float rsv_b_narrow_bandwidth = 0.1666667
    Float rsv_b_proportion_wide = 0.0
  }
  output {
    # standardized organism flag
    String standardized_organism = select_first([sc2_org_name, mpox_org_name, flu_org_name, rsv_a_org_name, rsv_b_org_name, organism])
    # reference genome and sequencing information
    File reference = select_first([reference_genome, sc2_reference_genome, mpox_reference_genome, h1n1_ha_reference, h3n2_ha_reference, vic_ha_reference, yam_ha_reference, h5n1_ha_reference, h1n1_na_reference, h3n2_na_reference, vic_na_reference, yam_na_reference, rsv_a_reference_genome, rsv_b_reference_genome, "gs://theiagen-public-resources-rp/empty_files/empty.fasta"])
    # augur options
    Int augur_min_num_unambig = select_first([min_num_unambig, sc2_min_num_unambig, mpox_min_num_unambig, flu_min_num_unambig, rsv_a_min_num_unambig, rsv_b_min_num_unambig, 0])
    File augur_clades_tsv = select_first([clades_tsv, sc2_clades_tsv, h1n1_ha_clades_tsv, h3n2_ha_clades_tsv, vic_ha_clades_tsv, yam_ha_clades_tsv, rsv_a_clades_tsv, rsv_b_clades_tsv, mpox_clades_tsv, "gs://theiagen-public-resources-rp/empty_files/minimal-clades.tsv"])
    File reference_gbk = select_first([reference_genbank, sc2_reference_genbank, h1n1_ha_reference_gbk, h3n2_ha_reference_gbk, vic_ha_reference_gbk, yam_ha_reference_gbk, h5n1_ha_reference_gbk, h1n1_na_reference_gbk, h3n2_na_reference_gbk, vic_na_reference_gbk, yam_na_reference_gbk, rsv_a_reference_gbk, rsv_b_reference_gbk, mpox_reference_gbk, "gs://theiagen-public-resources-rp/empty_files/empty.gbk"])
    File augur_auspice_config = select_first([auspice_config, sc2_auspice_config, h1n1_ha_auspice_config, h3n2_ha_auspice_config, vic_ha_auspice_config, yam_ha_auspice_config, h5n1_ha_auspice_config, h1n1_na_auspice_config, h3n2_na_auspice_config, vic_na_auspice_config, yam_na_auspice_config, rsv_a_auspice_config, rsv_b_auspice_config, mpox_auspice_config, "gs://theiagen-public-resources-rp/empty_files/minimal-auspice-config.json"])
    Float augur_min_date = select_first([min_date, sc2_min_date, flu_min_date, rsv_a_min_date, rsv_b_min_date, mpox_min_date, 0.0])
    Int augur_pivot_interval = select_first([pivot_interval, sc2_pivot_interval, flu_pivot_interval, mpox_pivot_interval, rsv_a_pivot_interval,rsv_b_pivot_interval, 0])
    Float augur_narrow_bandwidth = select_first([narrow_bandwidth, sc2_narrow_bandwidth, flu_narrow_bandwidth, mpox_narrow_bandwidth, rsv_a_narrow_bandwidth, rsv_b_narrow_bandwidth, 0.0])
    Float augur_proportion_wide = select_first([proportion_wide, sc2_proportion_wide, flu_proportion_wide,rsv_a_proportion_wide,rsv_b_proportion_wide,mpox_proportion_wide, 0.0])
  }
}
