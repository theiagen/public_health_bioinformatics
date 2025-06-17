# Augur

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Phylogenetic Construction](../../workflows_overview/workflows_type.md/#phylogenetic-construction) | [Viral](../../workflows_overview/workflows_kingdom.md/#viral) | PHB v3.0.0 | Yes | Sample-level, Set-level |

## Augur Workflows

Genomic Epidemiology is an important approach in the effort to understand and mitigate against disease transmission. An often-critical step in viral genomic epidemiology is the generation of phylogenetic trees to explore the genetic relationship between viruses on a local, regional, national or global scale. The Augur workflows, currently only targeted for viral pathogens, facilitate this process by generating files for the visualization of phylogenetic trees with accompanying metadata.

Two workflows are offered: **Augur_Prep_PHB** and **Augur_PHB**. These must be run sequentially, respectively, to first prepare each individual sample for running Augur, and secondly to run Augur itself on the set of samples, generating the phylogenetic tree files with accompanying metadata. The outputs from these workflows can be visualized in [Auspice](https://docs.nextstrain.org/projects/auspice/en/latest/) and [UShER](https://github.com/yatisht/usher).

!!! dna "**Helpful resources for epidemiological interpretation**"

    - [introduction to Nextstrain](https://www.cdc.gov/advanced-molecular-detection/php/training/module-3-1.html) (which includes Auspice)
    - guide to Nextstrain [interactive trees](https://www.cdc.gov/advanced-molecular-detection/php/training/module-3-4.html)
    - an [introduction to UShER](https://www.cdc.gov/advanced-molecular-detection/php/training/module-3-3.html)
    - a video about [how to read trees](https://www.cdc.gov/advanced-molecular-detection/php/training/module-1-3.html) if this is new to you
    - documentation on [how to identify SARS-CoV-2 recombinants](https://github.com/pha4ge/pipeline-resources/blob/main/docs/sc2-recombinants.md)

### Augur_Prep_PHB

The Augur_Prep_PHB workflow was written to prepare individual sample assemblies and their metadata for running the Augur_PHB analysis.

#### Augur_Prep Inputs

The Augur_Prep_PHB workflow takes assembly FASTA files and associated metadata formatted in a data table. FASTA files may be generated with one of the TheiaCoV Characterization workflows and should adhere to quality control guidelines, (e.g. [QC guidelines produced by PHA4GE](https://github.com/pha4ge/pipeline-resources/blob/main/docs/qc-solutions.md)). The metadata can be uploaded to Terra as TSV file, formatted as in [this example](https://docs.google.com/spreadsheets/d/1PF1u3R-ZGm53UiVsTlIcpg9Qk2dUJgtx/edit#gid=253517867).

This workflow runs on the sample level.

<div class="searchable-table" markdown="1">

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| augur_prep | **assembly** | File | Assembly/consensus file (single FASTA file per sample) |  | Required |
| augur_prep | **collection_date** | String | Collection date of the sample |  | Optional |
| augur_prep | **continent** | String | Continent where sample was collected |  | Optional |
| augur_prep | **country** | String | Country where sample was collected |  | Optional |
| augur_prep | **county** | String | County (or smaller locality) where sample was collected |  | Optional |
| augur_prep | **nextclade_clade** | String | The Nextclade clade of the sample |  | Optional |
| augur_prep | **pango_lineage** | String | The Pangolin lineage of the sample |  | Optional |
| augur_prep | state | **String** | State (or province) where sample was collected |  | Optional |
| prep_augur_metadata | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional |
| prep_augur_metadata | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 10 | Optional |
| prep_augur_metadata | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1 | Optional |
| prep_augur_metadata | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 3 | Optional |
| prep_augur_metadata | **organism** | String | The organism to be analyzed in Augur; options: "sars-cov-2", "flu", "MPXV", "rsv-a", "rsv-b" | sars-cov-2 | Optional |
| version_capture | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional |

</div>

#### Augur_Prep Outputs

| **Variable** | **Type** | **Description** |
|---|---|---|
| augur_metadata | File | TSV file of the metadata provided as input to the workflow in the proper format for Augur analysis |
| augur_prep_phb_analysis_date | String | Date of analysis |
| augur_prep_phb_version | String | Version of the Public Health Bioinformatics (PHB) repository used |

### Augur_PHB

!!! info "Helpful Hint"
    You may have to generate phylogenies multiple times, running the Augur_PHB workflow, assessing results, and amending inputs to generate a final tree with sufficient diversity and high-quality data of interest.

The Augur_PHB workflow takes a **set** of assembly/consensus files (FASTA format) and sample metadata files (TSV format) that have been reformatted using Augur_Prep_PHB and runs Augur to generate the phylogenetic tree files with accompanying metadata. Additionally, the workflow infers pairwise SNP distances.

#### Augur Inputs

The Augur_PHB workflow takes in a ***set*** of SARS-CoV-2 (or any other viral pathogen) FASTA and metadata files. If running the workflow via Terra, individual samples will need to be added to a set before running the workflow. Input FASTAs should meet QA metrics. Sets of FASTAs with highly discordant quality metrics may result in the inaccurate inference of genetic relatedness. There **must** be some sequence diversity among the set of input assemblies. If insufficient diversity is present, it may be necessary to add a more divergent sequence to the set.

!!! dna "Optional Inputs"
    There are **many** optional user inputs. For SARS-CoV-2, Flu, rsv-a, rsv-b, and mpxv, default values that mimic the NextStrain builds have been preselected. To use these defaults, you must write either `"sars-cov-2"`,`"flu"`, `"rsv-a"`, `"rsv-b"`, or `"mpxv"` for the `organism` variable.

    For Flu - it is **required** to set `flu_segment` to either `"HA"` or `"NA"` & `flu_subtype` to either `"H1N1"` or `"H3N2"` or `"Victoria"` or `"Yamagata"` or `"H5N1"` (`"H5N1"` will only work with `"HA"`) depending on your set of samples.

???+ toggle "A Note on Optional Inputs"
    ??? toggle "Default values for SARS-CoV-2"
        - min_num_unambig = 27000
        - clades_tsv = [defaults/clades.tsv](https://github.com/nextstrain/ncov/tree/23d1243127e8838a61b7e5c1a72bc419bf8c5a0d/defaults/clades.tsv)
        - lat_longs_tsv = [defaults/lat_longs.tsv](https://github.com/nextstrain/ncov/blob/23d1243127e8838a61b7e5c1a72bc419bf8c5a0d/defaults/lat_longs.tsv)
        - reference_fasta = [defaults/reference_seq.fasta](https://github.com/nextstrain/ncov/blob/23d1243127e8838a61b7e5c1a72bc419bf8c5a0d/defaults/reference_seq.fasta)
        - reference_genbank = [defaults/reference_seq.gb](https://github.com/nextstrain/ncov/blob/23d1243127e8838a61b7e5c1a72bc419bf8c5a0d/defaults/reference_seq.gb)
        - auspice_config = [defaults/auspice_config.json](https://github.com/nextstrain/ncov/blob/23d1243127e8838a61b7e5c1a72bc419bf8c5a0d/defaults/auspice_config.json)
        - min_date = 2020.0
        - pivot_interval = 1
        - pivot_interval_units = "weeks"
        - narrow_bandwidth = 0.05
        - proportion_wide = 0.0

    ??? toggle "Default values for Flu"
        - lat_longs_tsv = `"gs://theiagen-public-files-rp/terra/flu-references/lat_longs.tsv"`
        - min_num_unambig = 900
        - min_date = 2020.0
        - pivot_interval = 1
        - narrow_bandwidth = 0.1666667
        - proportion_wide = 0.0
        ??? toggle "H1N1"
            - auspice_config = `"gs://theiagen-public-files-rp/terra/flu-references/auspice_config_h1n1pdm.json"`
            - HA
                - reference_fasta = `"gs://theiagen-public-files-rp/terra/flu-references/reference_h1n1pdm_ha.gb"`
                - clades_tsv = `"gs://theiagen-public-files-rp/terra/flu-references/clades_h1n1pdm_ha.tsv"`
            - NA
                - reference_fasta = `"gs://theiagen-public-files-rp/terra/flu-references/reference_h1n1pdm_na.gb"`
        ??? toggle "H3N2"
            - auspice_config = `"gs://theiagen-public-files-rp/terra/flu-references/auspice_config_h3n2.json"`
            - HA
                - reference_fasta = `"gs://theiagen-public-files-rp/terra/flu-references/reference_h3n2_ha.gb"`
                - clades_tsv = `"gs://theiagen-public-files-rp/terra/flu-references/clades_h3n2_ha.tsv"`
            - NA
                - reference_fasta = `"gs://theiagen-public-files-rp/terra/flu-references/reference_h3n2_na.gb"`
        ??? toggle "Victoria"
            - auspice_config = `"gs://theiagen-public-files-rp/terra/flu-references/auspice_config_vic.json"`
            - HA
                - reference_fasta = `"gs://theiagen-public-files-rp/terra/flu-references/reference_vic_ha.gb"`
                - clades_tsv = `"gs://theiagen-public-files-rp/terra/flu-references/clades_vic_ha.tsv"`
            - NA
                - reference_fasta = `"gs://theiagen-public-files-rp/terra/flu-references/reference_vic_na.gb"`
        ??? toggle "Yamagata"
            - auspice_config = `"gs://theiagen-public-files-rp/terra/flu-references/auspice_config_yam.json"`
            - HA
                - reference_fasta = `"gs://theiagen-public-files-rp/terra/flu-references/reference_yam_ha.gb"`
                - clades_tsv = `"gs://theiagen-public-files-rp/terra/flu-references/clades_yam_ha.tsv"`
            - NA
                - reference_fasta = `"gs://theiagen-public-files-rp/terra/flu-references/reference_yam_na.gb"`
        ??? toggle "H5N1"
            - auspice_config = `"gs://theiagen-public-files-rp/terra/flu-references/auspice_config_h5n1.json"`
            - HA
                - reference_fasta = `"gs://theiagen-public-files-rp/terra/flu-references/reference_h5n1_ha.gb"`
                - clades_tsv = `"gs://theiagen-public-files-rp/terra/flu-references/h5nx-clades.tsv"`

    ??? toggle "Default values for MPXV"
        - min_num_unambig = 150000
        - clades_tsv = `"gs://theiagen-public-files-rp/terra/augur-mpox-references/mpox_clades.tsv"`
        - lat_longs_tsv = `"gs://theiagen-public-files-rp/terra/flu-references/lat_longs.tsv"`
        - reference_fasta = `"gs://theiagen-public-files-rp/terra/augur-mpox-references/NC_063383.1.reference.fasta"`
        - reference_genbank = `"gs://theiagen-public-files-rp/terra/augur-mpox-references/NC_063383.1_reference.gb"`
        - auspice_config = `"gs://theiagen-public-files-rp/terra/augur-mpox-references/mpox_auspice_config_mpxv.json"`
        - min_date = 2020.0
        - pivot_interval = 1
        - narrow_bandwidth = 0.1666667
        - proportion_wide = 0.0

    ??? toggle "Default values for RSV-A"
        - min_num_unambig = 10850
        - clades_tsv = `"gs://theiagen-public-files-rp/terra/rsv_references/rsv_a_clades.tsv"`
        - lat_longs_tsv = `"gs://theiagen-public-files-rp/terra/flu-references/lat_longs.tsv"`
        - reference_fasta = `"gs://theiagen-public-files-rp/terra/rsv_references/reference_rsv_a.fasta"`
        - reference_genbank = `""gs://theiagen-public-files-rp/terra/rsv_references/reference_rsv_a.gb"`
        - auspice_config = `""gs://theiagen-public-files-rp/terra/rsv_references/rsv_auspice_config.json"`
        - min_date = 2020.0
        - pivot_interval = 1
        - narrow_bandwidth = 0.1666667
        - proportion_wide = 0.0

    ??? toggle "Default values for RSV-B"
        - min_num_unambig = 10850
        - clades_tsv = `"gs://theiagen-public-files-rp/terra/rsv_references/rsv_b_clades.tsv"`
        - lat_longs_tsv = `"gs://theiagen-public-files-rp/terra/flu-references/lat_longs.tsv"`
        - reference_fasta = `"gs://theiagen-public-files-rp/terra/rsv_references/reference_rsv_b.fasta"`
        - reference_genbank = `""gs://theiagen-public-files-rp/terra/rsv_references/reference_rsv_b.gb"`
        - auspice_config = `""gs://theiagen-public-files-rp/terra/rsv_references/rsv_auspice_config.json"`
        - min_date = 2020.0
        - pivot_interval = 1
        - narrow_bandwidth = 0.1666667
        - proportion_wide = 0.0

    For more information regarding these optional inputs, please view [Nextrain's detailed documentation on Augur](https://docs.nextstrain.org/projects/augur/en/stable/usage/usage.html)

    !!! info "What's required or not?"
        For organisms _other_ than SARS-CoV-2 or Flu, the required variables have both the "required" and "optional" tags.

This workflow runs on the set level. Please note that for every task, runtime parameters are modifiable (cpu, disk_size, docker, and memory); most of these values have been excluded from the table below for convenience.

<div class="searchable-table" markdown="1"> 

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| augur | **assembly_fastas** | Array[File] | An array of the assembly files to use; use either the HA or NA segment for flu samples |  | Required |
| augur | **build_name** | String | Name to give to the Augur build |  | Required |
| augur | **auspice_config** | File | Auspice config file for customizing visualizations; takes priority over the other customization values available for augur_export | Defaults are organism-specific. Please find default values for all organisms (and for Flu - their respective genome segments and subtypes) here: <https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_organism_parameters.wdl>. For an organism without set defaults, a minimal auspice config file is provided to prevent workflow failure, "gs://theiagen-public-files-rp/terra/augur-defaults/minimal-auspice-config.json", but will not be as useful as an organism specific config file. | Optional |
| augur | **clades_tsv** | File | TSV file containing clade mutation positions in four columns | Defaults are organism-specific. Please find default values for all organisms (and for Flu - their respective genome segments and subtypes) here: <https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_organism_parameters.wdl>. For an organism without set defaults, an empty clades file is provided to prevent workflow failure, "gs://theiagen-public-files-rp/terra/augur-defaults/minimal-clades.tsv", but will not be as useful as an organism specific clades file. | Optional, Required |
| augur | **distance_tree_only** | Boolean | Create only a distance tree (skips all Augur steps after augur_tree) | TRUE | Optional |
| augur | **flu_segment** | String | Required if organism = "flu". The name of the segment to be analyzed; options: "HA" or "NA" | "HA" (only used if organism = "flu") | Optional, Required |
| augur | **flu_subtype** | String | Required if organism = "flu". The subtype of the flu samples being analyzed; options: "H1N1", "H3N2", "Victoria", "Yamagata", "H5N1" |  | Optional, Required |
| augur | **lat_longs_tsv** | File | Tab-delimited file of geographic location names with corresponding latitude and longitude values | Defaults are organism-specific. Please find default values for all organisms (and for Flu - their respective genome segments and subtypes) here: <https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_organism_parameters.wdl>. For an organism without set defaults, a minimal lat-long file is provided to prevent workflow failure, "gs://theiagen-public-files-rp/terra/augur-defaults/minimal-lat-longs.tsv", but will not be as useful as a detailed lat-longs file covering all the locations for the samples to be visualized. | Optional |
| augur | **min_date** | Float | Minimum date to begin filtering or frequencies calculations | Defaults are organism-specific. Please find default values for all organisms (and for Flu - their respective genome segments and subtypes) here: <https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_organism_parameters.wdl>. For an organism without set defaults, the default value is 0.0 | Optional |
| augur | **min_num_unambig** | Int | Minimum number of called bases in genome to pass prefilter | Defaults are organism-specific. Please find default values for all organisms (and for Flu - their respective genome segments and subtypes) here: <https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_organism_parameters.wdl>. For an organism without set defaults, the default value is 0 | Optional |
| augur | **organism** | String | Organism used to preselect default values; options: "sars-cov-2", "flu", "mpxv", "rsv-a", "rsv-b" | sars-cov-2 | Optional |
| augur | **reference_fasta** | File | The reference FASTA file used to align the genomes and build the trees | Defaults are organism-specific. Please find default values for all organisms (and for Flu - their respective genome segments and subtypes) here: <https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_organism_parameters.wdl>. For an organism without set defaults, a reference fasta file must be provided otherwise the workflow fails. | Optional, Required |
| augur | **reference_genbank** | File | The GenBank .gb file for the same reference genome used for the reference_fasta | Defaults are organism-specific. Please find default values for all organisms (and for Flu - their respective genome segments and subtypes) here: <https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_organism_parameters.wdl>. For an organism without set defaults, a reference genbank file must be provided otherwise the workflow fails. | Optional, Required |
| augur | **sample_metadata_tsvs** | Array[File] | An array of the metadata files produced in Augur_Prep_PHB |  | Optional |
| augur | **build_name_updated** | String | Internal component, do not modify. Used for replacing spaces with underscores _ |  | Do Not Modify |
| augur_align | **fill_gaps** | Boolean | If true, gaps represent missing data rather than true indels and so are replaced by N after aligning. | FALSE | Optional |
| augur_ancestral | **infer_ambiguous** | Boolean | If true, infer nucleotides and ambiguous sites and replace with most likely | FALSE | Optional |
| augur_ancestral | **inference** | String | Calculate joint or marginal maximum likelihood ancestral sequence states; options: "joint", "marginal" | joint | Optional |
| augur_ancestral | **keep_ambiguous** | Boolean | If true, do not infer nucleotides at ambiguous (N) sides | FALSE | Optional |
| augur_ancestral | **keep_overhangs** | Boolean | If true, do not infer nucleotides for gaps on either side of the alignment | FALSE | Optional |
| augur_export | **colors_tsv** | File | Custom color definitions, one per line in TSV format with the following fields: TRAIT_TYPE TRAIT_VALUE HEX_CODE |  | Optional |
| augur_export | **description_md** | File | Markdown file with description of build and/or acknowledgements |  | Optional |
| augur_export | **include_root_sequence** | Boolean | Export an additional JSON containing the root sequence used to identify mutations | FALSE | Optional |
| augur_export | **title** | String | Title to be displayed by Auspice |  | Optional |
| augur_refine | **branch_length_inference** | String | Branch length mode of timetree to use; options: "auto", "joint", "marginal", "input" | auto | Optional |
| augur_refine | **clock_filter_iqd** | Int | Remove tips that deviate more than n_iqd interquartile ranges from the root-to-tip vs time regression | 4 | Optional |
| augur_refine | **clock_rate** | Float | Fixed clock rate to use for time tree calculations |  | Optional |
| augur_refine | **clock_std_dev** | Float | Standard deviation of the fixed clock_rate estimate |  | Optional |
| augur_refine | **coalescent** | String | Coalescent time scale in units of inverse clock rate (float), optimize as scalar ("opt") or skyline ("skyline") |  | Optional |
| augur_refine | **covariance** | Boolean | If true, account for covariation when estimating rates and/or rerooting | TRUE | Optional |
| augur_refine | **date_confidence** | Boolean | If true, calculate confidence intervals for node dates | TRUE | Optional |
| augur_refine | **date_inference** | String | Assign internal nodes to their marginally most likely dates; options: "joint", "marginal" | marginal | Optional |
| augur_refine | **divergence_units** | String | Units in which sequence divergences is exported; options: "mutations" or "mutations-per-site" | mutations | Optional |
| augur_refine | **gen_per_year** | Int | Number of generations per year | 50 | Optional |
| augur_refine | **keep_polytomies** | Boolean | If true, don't attempt to resolve polytomies | FALSE | Optional |
| augur_refine | **keep_root** | Boolean | If true, do not reroot the tree; use it as-is (overrides anything specified by root) | TRUE | Optional |
| augur_refine | **precision** | String | Precision used to determine the number of grid points; options: 0 (rough) to 3 (ultra fine) | auto | Optional |
| augur_refine | **root** | String | Rooting mechanism; options: "best", "least-squares", "min_dev", "oldest", etc. |  | Optional |
| augur_translate | **genes** | File | A file containing a list of genes to translate (from nucleotides to amino acids) |  | Optional |
| augur_tree | **exclude_sites** | File | File of one-based sites to exclude for raw tree building (BED format in .bed files, DRM format in tab-delimited files, or one position per line) |  | Optional |
| augur_tree | **method** | String | Which method to use to build the tree; options: "fasttree", "raxml", "iqtree" | iqtree | Optional |
| augur_tree | **override_default_args** | Boolean | If true, override default tree builder arguments instead of augmenting them | FALSE | Optional |
| augur_tree | **substitution_model** | String | The substitution model to use; only available for iqtree. Specify "auto" to run ModelTest; model options can be found [here](http://www.iqtree.org/doc/Substitution-Models) | GTR | Optional |
| augur_tree | **tree_builder_args** | String | Additional tree builder arguments either augmenting or overriding the default arguments. FastTree defaults: "-nt -nosupport". RAxML defaults: "-f d -m GTRCAT -c 25 -p 235813". IQ-TREE defaults: "-ninit 2 -n 2 -me 0.05 -nt AUTO -redo" |  | Optional |
| sc2_defaults | **nextstrain_ncov_repo_commit** | String | The version of the <https://github.com/nextstrain/ncov/> from which to draw default values for SARS-CoV-2. | `23d1243127e8838a61b7e5c1a72bc419bf8c5a0d` | Optional |
| organism_parameters | **gene_locations_bed_file** | File | Use to provide locations of interest where average coverage will be calculated | Defaults are organism-specific. Please find default values for some organisms here: <https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_organism_parameters.wdl>. For an organism without set defaults, an empty file is provided, "gs://theiagen-public-files/terra/theiacov-files/empty.bed", but will not be as useful as an organism specific gene locations bed file. | Optional |
| organism_parameters | **genome_length_input** | Int | Use to specify the expected genome length; provided by default for all supported organisms | Defaults are organism-specific. Please find default values for all organisms (and for Flu - their respective genome segments and subtypes) here: <https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_organism_parameters.wdl>. For an organism without set defaults, the genome length input must be provided otherwise the workflow fails. | Optional, Required |
| organism_parameters | **hiv_primer_version** | String | The version of HIV primers used. Options are <https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_organism_parameters.wdl#L156> and <https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_organism_parameters.wdl#L164>. This input is ignored if provided for TheiaCoV_Illumina_SE and TheiaCoV_ClearLabs | v1 | Optional |
| organism_parameters | **kraken_target_organism_input** | String | The organism whose abundance the user wants to check in their reads. This should be a proper taxonomic name recognized by the Kraken database. | Defaults are organism-specific. Please find default values for all organisms here: <https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_organism_parameters.wdl>. For an organism without set defaults, the default is "". | Optional |
| organism_parameters | **nextclade_dataset_name_input** | String | NextClade organism dataset name | Defaults are organism-specific. Please find default values for all organisms (and for Flu - their respective genome segments and subtypes) here: <https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_organism_parameters.wdl>. For an organism without set defaults, the default is "NA". | Optional |
| organism_parameters | **nextclade_dataset_tag_input** | String | NextClade organism dataset tag | Defaults are organism-specific. Please find default values for all organisms (and for Flu - their respective genome segments and subtypes) here: <https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_organism_parameters.wdl>. For an organism without set defaults, the default is "NA". | Optional |
| organism_parameters | **pangolin_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/pangolin:4.3.1-pdata-1.26 | Optional |
| organism_parameters | **primer_bed_file** | File | The bed file containing the primers used when sequencing was performed | Defaults are organism-specific. Please find default values for all organisms here: <https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_organism_parameters.wdl>. For an organism without set defaults, an empty primer bed file is provided, "gs://theiagen-public-files/terra/theiacov-files/empty.bed", but will not be as useful as an organism specific primer bed file. | Optional |
| organism_parameters | **reference_gff_file** | File | Reference GFF file for the organism being analyzed | Defaults are organism-specific. Please find default values for all organisms here: <https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_organism_parameters.wdl>. For an organism without set defaults, an empty gff file is provided, "gs://theiagen-public-files/terra/theiacov-files/empty.gff3", but will not be as useful as an organism specific gff file. | Optional |
| organism_parameters | **vadr_max_length** | Int | Maximum length for the `fasta-trim-terminal-ambigs.pl` VADR script | Defaults are organism-specific. Please find default values for all organisms here: <https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_organism_parameters.wdl>. For an organism without set defaults, the default is 0. | Optional |
| organism_parameters | **vadr_mem** | Int | Memory, in GB, allocated to this task | 32 (RSV-A and RSV-B) and 8 (all other TheiaCoV organisms) |  |
| organism_parameters | **vadr_options** | String | Options for the `v-annotate.pl` VADR script | Defaults are organism-specific. Please find default values for all organisms here: <https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_organism_parameters.wdl>. For an organism without set defaults, the default is "NA". | Optional |
| organism_parameters | **vadr_skip_length** | Int | Minimum assembly length (unambiguous) to run VADR | Defaults are organism-specific. Please find default values for all organisms here: <https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_organism_parameters.wdl>. For an organism without set defaults, the default is 0. | Optional |
| mutation_context | **cpu** | Int | CPUs requested for the mutation_context task that is specific to Mpox. | 1 | Optional |
| mutation_context | **disk_size** | Int | Disk size in GB requested for the mutation_context task that is specific to Mpox. | 50 | Optional |
| mutation_context | **docker** | String | Docker image used for the mutation_context task that is specific to Mpox. Do not modify. | us-docker.pkg.dev/general-theiagen/theiagen/nextstrain-mpox-mutation-context:2024-06-27 | Do Not Modify, Optional |
| mutation_context | **memory** | Int | Memory size in GB requested for the mutation_context task that is specific to Mpox. | 4 | Optional |

</div>

??? task "Workflow Tasks"
    ##### Augur Workflow Tasks {#augur-tasks}

    The Augur_PHB workflow uses the inputs to generate a phylogenetic tree in JSON format that is compatible with phylogenetic tree visualization software.
    
    In Augur_PHB, the tasks below are called. For the Augur subcommands, please view the [Nextstrain Augur documentation](https://docs.nextstrain.org/projects/augur/en/stable/usage/usage.html) for more details and explanations.

    1. `cat_files` - concatenate all of the input fasta files together
    2. `sc2_defaults` - if organism is SARS-CoV-2, establish default parameters
    3. `flu_defaults` - if organism is Flu, establish default parameters
    4. `filter_sequences_by_length` - remove any sequences that do not meet the quality threshold set by `min_num_unambig`
    5. `tsv_join` - merge the metadata files
    6. `fasta_to_ids` - extract a list of remaining sequences so we know which ones were dropped
    7. `augur_align` - perform MAFFT alignment on the sequences
    8. `augur_tree` - create a distance tree
    9. `augur_refine` - create a timetree
    10. `augur_ancestral` - infer ancestral sequences
    11. `augur_translate` - translate gene regions from nucleotides to amino acids
    12. `mutation_context` - if organism is MPXV, calculates the mutation fraction of G->A or C->T changes
    13. `augur_clades` - if clade information is provided, assign clades to nodes based on amino-acid or nucleotide signatures
    14. `augur_export` - export all the results in a JSON file suitable for Auspice visualization
    15. `snp_dists` - create a SNP matrix from the alignment
    16. `reorder_matrix` - reorder the SNP matrix to match the distance tree

#### Augur Outputs

!!! dna "Diversity dependent"
    Note that the node & branch coloring by clade or lineage assignment might be dependent on the diversity of your input dataset. This is because the clade assignment is done using the ancestrally reconstructed amino acid or nucleotide changes at the tree nodes rather than a direct sequence-to-reference mutation comparison. You may notice this happening when you get clade/lineage assignments from NextClade when running TheiaCoV workflows, but no clade/lineage assignment on the Augur Auspice tree.

    To get around this issue, you can upload the Augur output file `merged-metadata.tsv` to Auspice that includes the correct clade/lineage assignments to allow for coloring by Clade.

!!! dna "Flu clade assignments"
    Note that for flu, the clade assignment is usually mostly done for the more recent seasonal influenza viruses. Older strains may get an "unassigned" designation for clades. Therefore, it is important to counter check with the NextClade results from TheiaCoV if the lack of clade assignment is due to analyzing older sequences or sequence related.

The `auspice_input_json` is intended to be uploaded to [Auspice](https://auspice.us/) to view the phylogenetic tree. This provides a visualization of the genetic relationships between your set of samples. The `metadata_merged` output can also be uploaded to add context to the phylogenetic visualization. The `combined_assemblies` output can be uploaded to [UShER](https://genome.ucsc.edu/cgi-bin/hgPhyloPlace) to view the samples on a global tree of representative sequences from the public repositories.

The Nextstrain team hosts documentation surrounding the Augur workflow → Auspice visualization here, which details the various components of the Auspice interface: [How data is exported by Augur for visualisation in Auspice](https://docs.nextstrain.org/en/latest/learn/augur-to-auspice.html).

| **Variable** | **Type** | **Description** |
| --- | --- | --- |
| aligned_fastas | File | A FASTA file of the aligned genomes |
| augur_fasttree_version | String | The fasttree version used, blank if other tree method used |
| augur_iqtree_model_used | String | The iqtree model used during augur tree, blank if iqtree not used |
| augur_iqtree_version | String | The iqtree version used during augur tree (defualt), blank if other tree method used |
| augur_mafft_version | String | The mafft version used in augur align |
| augur_phb_analysis_date | String | The date the analysis was run |
| augur_phb_version | String | The version of the Public Health Bioinformatics (PHB) repository used |
| augur_raxml_version | String | The version of raxml used during augur tree, blank if other tree method used |
| augur_version | String | Version of Augur used |
| auspice_input_json | File | JSON file used as input to Auspice |
| combined_assemblies | File | Concatenated FASTA file containing all samples |
| distance_tree | File | The distance tree created in Newick (.nwk) format |
| keep_list | File | A list of samples included in the phylogenetic tree |
| metadata_merged | File | Tab-delimited text file of the merged augur_metadata input files from all samples |
| snp_matrix | File | The SNP distance matrix for all samples used in the phylogenetic tree |
| time_tree | File | The time tree created in Newick (.nwk) format |
| traits_json | File | A JSON file containing sample traits |

#### Mpox-specific Auspice Output JSON

If you are building a tree for Mpox samples and set the optional input parameter `organism` to `"mpox"` , an additional step will be carried out in the Augur_PHB workflow. This additional step will calculate the mutation fraction of G→A or C→T changes. These mutations have been shown to be a characteristic of APOBEC3-type editing, which indicate adaptation of the virus to circulation among humans as was observed with the 2022 clade IIb outbreak, and more recently (2024) with the clade Ib outbreak in South Kivu, Democratic Republic of the Congo.

When visualizing the output `auspice_input_json` file, there will be 2 new choices in the drop-down menu for "Color By":

- G→A or C→T fraction
- NGA/TCN context of G→A or C→T mutations.

An example Mpox tree with these "Color By" options can be viewed here: <https://nextstrain.org/mpox/clade-IIb?c=GA_CT_fraction>

### References

When publishing work using the Augur_PHB workflow, please reference the following:

> Nextstrain: Hadfield J, Megill C, Bell SM, Huddleston J, Potter B, Callender C, Sagulenko P, Bedford T, Neher RA. Nextstrain: real-time tracking of pathogen evolution. Bioinformatics. 2018 Dec 1;34(23):4121-3.

When publishing work using inferences from UShER, please reference:

> UShER: Turakhia Y, Thornlow B, Hinrichs AS, De Maio N, Gozashti L, Lanfear R, Haussler D, Corbett-Detig R. Ultrafast Sample placement on Existing tRees (UShER) enables real-time phylogenetics for the SARS-CoV-2 pandemic. Nature Genetics. 2021 Jun;53(6):809-16.
