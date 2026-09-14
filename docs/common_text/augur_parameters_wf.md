---
title: Workflow Fragment `augur_parameters`
fragment: true
---
**The `augur_parameters` sub-workflow is the first step in the Augur workflow**. This step automatically sets the reference files and Augur parameters to the appropriate value for the user-designated organism (`"sars-cov-2"` is the default organism).

!!! dna ""
    The following tables include the relevant organism-specific parameters; **all of these default values can be overwritten by providing a value for the "Overwrite Variable Name" field**.

    === "SARS-CoV-2"
        <div class="searchable-table" markdown="block">

        | **Overwrite Variable Name** | **Organism** | **Default Value** |
        |---|---|---|
        | min_num_unambig | sars-cov-2 | `27000` |
        | clades_tsv | sars-cov-2 | `"gs://theiagen-public-resources-rp/reference_data/viral/sars-cov-2/sc2_clades_20251008.tsv"` |
        | auspice_config | sars-cov-2 | `"gs://theiagen-public-resources-rp/reference_data/viral/sars-cov-2/sc2_auspice_config_20251030.json"` |
        | reference_fasta | sars-cov-2 | `"gs://theiagen-public-resources-rp/reference_data/viral/sars-cov-2/MN908947.fasta"` |
        | reference_genbank | sars-cov-2 | `"gs://theiagen-public-resources-rp/reference_data/viral/sars-cov-2/sc2_reference_seq_20251008.gb"` |

        </div>

    === "Mpox"
        <div class="searchable-table" markdown="block">

        | **Overwrite Variable Name** | **Organism** | **Default Value** |
        |---|---|---|
        | min_num_unambig | MPXV | `150000` |
        | clades_tsv | MPXV | `"gs://theiagen-public-resources-rp/reference_data/viral/mpox/mpox_clades.tsv"` |
        | auspice_config | MPXV | `"gs://theiagen-public-resources-rp/reference_data/viral/mpox/mpox_auspice_config_mpxv.json"` |
        | reference_fasta | MPXV | `"gs://theiagen-public-resources-rp/reference_data/viral/mpox/MPXV.MT903345.reference.fasta"` |
        | reference_genbank | MPXV | `"gs://theiagen-public-resources-rp/reference_data/viral/mpox/NC_063383.1_reference.gb"` |

        </div>

    === "Influenza"
        !!! dna "Defaults available for several subtypes"
            We have default Augur parameters available for the following subtypes:

            - H1N1 (HA and NA segments)
            - H3N2 (HA and NA segments)
            - Victoria (HA and NA segments)
            - Yamagata (HA and NA segments)
            - H5Nx (HA segment)

            A default `clades_tsv` is only available for the HA segment of H1N1, H3N2, Victoria, and Yamagata. You can find more details below.

        <div class="searchable-table" markdown="block">

        | **Overwrite Variable Name** | **Organism** | **Flu Segment** | **Flu Subtype** | **Default Value** |
        |---|---|---|---|---|
        | min_num_unambig | flu | all | all | `900` |
        | auspice_config | flu | ha | h1n1 | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/auspice_config_h1n1pdm.json"` |
        | clades_tsv | flu | ha | h1n1 | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/clades_h1n1pdm_ha.tsv"` |
        | reference_fasta | flu | ha | h1n1 | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h1n1pdm_ha.fasta"` |
        | reference_genbank | flu | ha | h1n1 | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h1n1pdm_ha.gb"` |
        | auspice_config | flu | ha | h3n2 | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/auspice_config_h3n2.json"` |
        | clades_tsv | flu | ha | h3n2 | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/clades_h3n2_ha.tsv"` |
        | reference_fasta | flu | ha | h3n2 | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h3n2_ha.fasta"` |
        | reference_genbank | flu | ha | h3n2 | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h3n2_ha.gb"` |
        | auspice_config | flu | ha | victoria | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/auspice_config_vic.json"` |
        | clades_tsv | flu | ha | victoria | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/clades_vic_ha.tsv"` |
        | reference_fasta | flu | ha | victoria | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_vic_ha.fasta"` |
        | reference_genbank | flu | ha | victoria | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_vic_ha.gb"` |
        | auspice_config | flu | ha | yamagata | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/auspice_config_yam_20251030.json"` |
        | clades_tsv | flu | ha | yamagata | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/clades_yam_ha.tsv"` |
        | reference_fasta | flu | ha | yamagata | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_yam_ha.fasta"` |
        | reference_genbank | flu | ha | yamagata | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_yam_ha.gb"` |
        | auspice_config | flu | ha | h5n1 | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/auspice_config_h5n1_20251030.json"` |
        | reference_fasta | flu | ha | h5n1 | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h5n1_ha.fasta"` |
        | reference_genbank | flu | ha | h5n1 | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h5n1_ha.gb"` |
        | auspice_config | flu | na | h1n1 | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/auspice_config_h1n1pdm_20251030.json"` |
        | reference_fasta | flu | na | h1n1 | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h1n1pdm_na.fasta"` |
        | reference_genbank | flu | na | h1n1 | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h1n1pdm_na.gb"` |
        | auspice_config | flu | na | h3n2 | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/auspice_config_h3n2_20251030.json"` |
        | reference_fasta | flu | na | h3n2 | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h3n2_na.fasta"` |
        | reference_genbank | flu | na | h3n2 | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h3n2_na.gb"` |
        | auspice_config | flu | na | victoria | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/auspice_config_vic_20251030.json"` |
        | reference_fasta | flu | na | victoria | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_vic_na.fasta"` |
        | reference_genbank | flu | na | victoria | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_vic_na.gb"` |
        | auspice_config | flu | na | yamagata | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/auspice_config_yam_20251030.json"` |
        | reference_fasta | flu | na | yamagata | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_yam_na.fasta"` |
        | reference_genbank | flu | na | yamagata | `"gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_yam_na.gb"` |

        </div>

        !!! tip "H5 subtypes"
            Any subtype matching `H5N*` (and `"H5"`) uses the H5N1 defaults.

    === "RSV-A"
        <div class="searchable-table" markdown="block">

        | **Overwrite Variable Name** | **Organism** | **Default Value** | **Notes** |
        |---|---|---|---|
        | min_num_unambig | rsv_a | `10850` | 70% of the 15500 bp RSV genome |
        | clades_tsv | rsv_a | `"gs://theiagen-public-resources-rp/reference_data/viral/rsv/rsv_a_clades.tsv"` | |
        | auspice_config | rsv_a | `"gs://theiagen-public-resources-rp/reference_data/viral/rsv/rsv_auspice_config_20251030.json"` | |
        | reference_fasta | rsv_a | `"gs://theiagen-public-resources-rp/reference_data/viral/rsv/reference_rsv_a.EPI_ISL_412866.fasta"` | |
        | reference_genbank | rsv_a | `"gs://theiagen-public-resources-rp/reference_data/viral/rsv/reference_rsv_a.gb"` | |

        </div>

    === "RSV-B"
        <div class="searchable-table" markdown="block">

        | **Overwrite Variable Name** | **Organism** | **Default Value** | **Notes** |
        |---|---|---|---|
        | min_num_unambig | rsv_b | `10850` | 70% of the 15500 bp RSV genome |
        | clades_tsv | rsv_b | `"gs://theiagen-public-resources-rp/reference_data/viral/rsv/rsv_b_clades.tsv"` | |
        | auspice_config | rsv_b | `"gs://theiagen-public-resources-rp/reference_data/viral/rsv/rsv_auspice_config_20251030.json"` | |
        | reference_fasta | rsv_b | `"gs://theiagen-public-resources-rp/reference_data/viral/rsv/reference_rsv_b.EPI_ISL_1653999.fasta"` | |
        | reference_genbank | rsv_b | `"gs://theiagen-public-resources-rp/reference_data/viral/rsv/reference_rsv_b.gb"` | |

        </div>

    !!! techdetails "Augur Parameters Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [wf_augur_parameters.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_augur_parameters.wdl) |
