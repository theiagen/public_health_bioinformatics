??? task "`organism_parameters`: Setting default values for specific organisms"

    Organism Parameters acquires and propagates default files and variable values for specific organisms.

    <!--if: virus -->

    ??? toggle "Default values for SARS-CoV-2"
        - min_num_unambig = 27000
        - clades_tsv = `"gs://theiagen-public-resources-rp/reference_data/viral/sars-cov-2/sc2_clades_20251008.tsv"`
        - lat_longs_tsv = `"gs://theiagen-public-resources-rp/reference_data/viral/sars-cov-2/sc2_lat_longs_20251008.tsv"`
        - reference_fasta = `"gs://theiagen-public-resources-rp/reference_data/viral/sars-cov-2/MN908947.fasta"`
        - reference_genbank = `"gs://theiagen-public-resources-rp/reference_data/viral/sars-cov-2/sc2_reference_seq_20251008.gb"`
        - auspice_config = `"gs://theiagen-public-resources-rp/reference_data/viral/sars-cov-2/sc2_auspice_config_20251008.json"`
        - min_date = 2020.0
        - pivot_interval = 1
        - pivot_interval_units = "weeks"
        - narrow_bandwidth = 0.05
        - proportion_wide = 0.0

    ??? toggle "Default values for Flu"
        - lat_longs_tsv = `"gs://theiagen-public-resources-rp/reference_data/viral/flu/lat_longs.tsv"`
        - min_num_unambig = 900
        - min_date = 2020.0
        - pivot_interval = 1
        - narrow_bandwidth = 0.1666667
        - proportion_wide = 0.0
        ??? toggle "H1N1"
            - auspice_config = `"gs://theiagen-public-resources-rp/reference_data/viral/flu/auspice_config_h1n1pdm.json"`
            - HA
                - reference_fasta = `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h1n1pdm_ha.gb"`
                - clades_tsv = `"gs://theiagen-public-resources-rp/reference_data/viral/flu/clades_h1n1pdm_ha.tsv"`
            - NA
                - reference_fasta = `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h1n1pdm_na.gb"`
        ??? toggle "H3N2"
            - auspice_config = `"gs://theiagen-public-resources-rp/reference_data/viral/flu/auspice_config_h3n2.json"`
            - HA
                - reference_fasta = `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h3n2_ha.gb"`
                - clades_tsv = `"gs://theiagen-public-resources-rp/reference_data/viral/flu/clades_h3n2_ha.tsv"`
            - NA
                - reference_fasta = `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h3n2_na.gb"`
        ??? toggle "Victoria"
            - auspice_config = `"gs://theiagen-public-resources-rp/reference_data/viral/flu/auspice_config_vic.json"`
            - HA
                - reference_fasta = `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_vic_ha.gb"`
                - clades_tsv = `"gs://theiagen-public-resources-rp/reference_data/viral/flu/clades_vic_ha.tsv"`
            - NA
                - reference_fasta = `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_vic_na.gb"`
        ??? toggle "Yamagata"
            - auspice_config = `"gs://theiagen-public-resources-rp/reference_data/viral/flu/auspice_config_yam.json"`
            - HA
                - reference_fasta = `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_yam_ha.gb"`
                - clades_tsv = `"gs://theiagen-public-resources-rp/reference_data/viral/flu/clades_yam_ha.tsv"`
            - NA
                - reference_fasta = `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_yam_na.gb"`
        ??? toggle "H5N1"
            - auspice_config = `"gs://theiagen-public-resources-rp/reference_data/viral/flu/auspice_config_h5n1.json"`
            - HA
                - reference_fasta = `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h5n1_ha.gb"`
                - clades_tsv = `"gs://theiagen-public-resources-rp/reference_data/viral/flu/h5nx-clades.tsv"`

    ??? toggle "Default values for MPXV"
        - min_num_unambig = 150000
        - clades_tsv = `"gs://theiagen-public-resources-rp/reference_data/viral/mpox/mpox_clades.tsv"`
        - lat_longs_tsv = `"gs://theiagen-public-resources-rp/reference_data/viral/flu/lat_longs.tsv"`
        - reference_fasta = `"gs://theiagen-public-resources-rp/reference_data/viral/mpox/NC_063383.1.reference.fasta"`
        - reference_genbank = `"gs://theiagen-public-resources-rp/reference_data/viral/mpox/NC_063383.1_reference.gb"`
        - auspice_config = `"gs://theiagen-public-resources-rp/reference_data/viral/mpox/mpox_auspice_config_mpxv.json"`
        - min_date = 2020.0
        - pivot_interval = 1
        - narrow_bandwidth = 0.1666667
        - proportion_wide = 0.0

    ??? toggle "Default values for RSV-A"
        - min_num_unambig = 10850
        - clades_tsv = `"gs://theiagen-public-resources-rp/reference_data/viral/rsv/rsv_a_clades.tsv"`
        - lat_longs_tsv = `"gs://theiagen-public-resources-rp/reference_data/viral/flu/lat_longs.tsv"`
        - reference_fasta = `"gs://theiagen-public-resources-rp/reference_data/viral/rsv/reference_rsv_a.EPI_ISL_412866.fasta"`
        - reference_genbank = `""gs://theiagen-public-resources-rp/reference_data/viral/rsv/reference_rsv_a.gb"`
        - auspice_config = `""gs://theiagen-public-resources-rp/reference_data/viral/rsv/rsv_auspice_config.json"`
        - min_date = 2020.0
        - pivot_interval = 1
        - narrow_bandwidth = 0.1666667
        - proportion_wide = 0.0

    ??? toggle "Default values for RSV-B"
        - min_num_unambig = 10850
        - clades_tsv = `"gs://theiagen-public-resources-rp/reference_data/viral/rsv/rsv_b_clades.tsv"`
        - lat_longs_tsv = `"gs://theiagen-public-resources-rp/reference_data/viral/flu/lat_longs.tsv"`
        - reference_fasta = `"gs://theiagen-public-resources-rp/reference_data/viral/rsv/reference_rsv_b.EPI_ISL_1653999.fasta"`
        - reference_genbank = `""gs://theiagen-public-resources-rp/reference_data/viral/rsv/reference_rsv_b.gb"`
        - auspice_config = `""gs://theiagen-public-resources-rp/reference_data/viral/rsv/rsv_auspice_config.json"`
        - min_date = 2020.0
        - pivot_interval = 1
        - narrow_bandwidth = 0.1666667
        - proportion_wide = 0.0

    <!-- endif -->

  
    !!! techdetails "Organism Parameters Technical Details"        

        |  | Links |
        | --- | --- |
        | Task | [wf_organism_parameters.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_organism_parameters.wdl) |
