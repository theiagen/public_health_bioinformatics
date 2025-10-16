# Augur

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Name": "[**Augur**](../workflows/phylogenetic_construction/augur.md)"}, columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level", "Dockstore"]) }}

## Augur Workflows

Genomic Epidemiology is an important approach to understand and mitigate disease transmission. A critical step in viral genomic epidemiology is generating phylogenetic trees to explore the genetic relationship between viruses on a local, regional, national, or global scale. To this end, the Augur workflows enable viral phylogenetic analysis by generating phylogenetic trees and incorporating metadata into a stunning visual via the Auspice format.

Two workflows are offered: **Augur_Prep_PHB** and **Augur_PHB**. These workflows must be run sequentially to generate a phylogenetic tree with metadata included, though **Augur_PHB** can generate a phylogenetic tree without metadata as a standalone workflow. The outputs from these workflows can be visualized in [Auspice](https://docs.nextstrain.org/projects/auspice/en/latest/) and [UShER](https://github.com/yatisht/usher).

!!! dna "**Helpful resources for epidemiological interpretation**"

    - [introduction to Nextstrain](https://www.cdc.gov/advanced-molecular-detection/php/training/module-3-1.html) (which includes Auspice)
    - guide to Nextstrain [interactive trees](https://www.cdc.gov/advanced-molecular-detection/php/training/module-3-4.html)
    - an [introduction to UShER](https://www.cdc.gov/advanced-molecular-detection/php/training/module-3-3.html)
    - a video about [how to read trees](https://www.cdc.gov/advanced-molecular-detection/php/training/module-1-3.html) if this is new to you
    - documentation on [how to identify SARS-CoV-2 recombinants](https://github.com/pha4ge/pipeline-resources/blob/main/docs/sc2-recombinants.md)

### Augur_Prep_PHB

The Augur_Prep_PHB workflow was written to prepare individual sample assemblies and their metadata for running the Augur_PHB analysis.

#### Augur_Prep Inputs

The Augur_Prep_PHB workflow takes assembly FASTA files and associated metadata formatted in a data table. FASTA files may be generated with one of the TheiaCoV/TheiaViral Characterization workflows and should adhere to quality control guidelines, (e.g. [QC guidelines produced by PHA4GE](https://github.com/pha4ge/pipeline-resources/blob/main/docs/qc-solutions.md)). The metadata can be uploaded to Terra as TSV file, formatted as in [this example](https://docs.google.com/spreadsheets/d/1PF1u3R-ZGm53UiVsTlIcpg9Qk2dUJgtx/edit#gid=253517867).

This workflow runs on the sample level.

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "Augur_Prep"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

#### Augur_Prep Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "Augur_Prep"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///

### Augur_PHB

!!! info "Helpful Hint"
    You may have to generate phylogenies multiple times, running the Augur_PHB workflow, assessing results, and amending inputs to generate a final tree with sufficient diversity and high-quality data of interest.

The Augur_PHB workflow takes a **set** of assembly/consensus files (FASTA format) and optional sample metadata files (TSV format) that have been reformatted using Augur_Prep_PHB. Augur_PHB runs Augur to generate the phylogenetic tree files with accompanying metadata. Additionally, the workflow infers pairwise SNP distances and can generate a time-calibrated phylogeny if collection date metadata was included in Augur_Prep_PHB.

#### Augur Inputs

The Augur_PHB workflow takes in a ***set*** of viral pathogen FASTA and metadata files. If running the workflow via Terra, individual samples will need to be added to a set before running the workflow. Input FASTAs should meet QC metrics. Sets of FASTAs with highly discordant quality metrics may result in the inaccurate inference of genetic relatedness. There **must** be some sequence diversity among the set of input assemblies. If insufficient diversity is present, it may be necessary to add a more divergent sequences to the set.

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

    For more information regarding these optional inputs, please view [Nextrain's detailed documentation on Augur](https://docs.nextstrain.org/projects/augur/en/stable/usage/usage.html)

    !!! info "What's required or not?"
        For organisms _other_ than SARS-CoV-2 or Flu, the required variables have both the "required" and "optional" tags.

This workflow runs on the set level. Please note that for every task, runtime parameters are modifiable (cpu, disk_size, docker, and memory); most of these values have been excluded from the table below for convenience.

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "Augur"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Workflow Tasks

The Augur_PHB workflow uses the inputs to generate a phylogenetic tree in Auspice JSON and Newick formats. 
    
In Augur_PHB, the tasks below are called. For the Augur subcommands, please view the [Nextstrain Augur documentation](https://docs.nextstrain.org/projects/augur/en/stable/usage/usage.html) for more details and explanations.

??? toggle "Versioning"

{{ include_md("common_text/versioning_task.md", indent=8) }}

??? toggle "Organism Parameters"

{{ include_md("common_text/organism_parameters_wf.md", indent=8) }}

??? toggle "Join Metadata TSVs"

{{ include_md("common_text/tsv_join_task.md", indent=8) }}

??? toggle "Concatenate Assemblies"

{{ include_md("common_text/cat_files_task.md", indent=8) }}

??? toggle "Filter Sequences by Length"

{{ include_md("common_text/filter_contigs_task.md", indent=8) }}

??? toggle "Augur Align"

{{ include_md("common_text/augur_align_task.md", indent=8) }}

??? toggle "Extract Sequence IDs"

{{ include_md("common_text/fasta_to_ids_task.md", indent=8) }}

??? toggle "Augur Tree"

{{ include_md("common_text/augur_tree_task.md", indent=8) }}

??? toggle "Calculate SNP Distances"

{{ include_md("common_text/snp_dists_task.md", indent=8) }}

??? toggle "Root and Reorder Matrix"

{{ include_md("common_text/reorder_matrix_task.md", indent=8) }}

??? toggle "Augur Refine"

{{ include_md("common_text/augur_refine_task.md", indent=8) }}

??? toggle "Augur Ancestral"

{{ include_md("common_text/augur_ancestral_task.md", indent=8) }}

??? toggle "Augur Translate"

{{ include_md("common_text/augur_translate_task.md", indent=8) }}

??? toggle "Add Mutation Context"

{{ include_md("common_text/mutation_context_task.md", indent=8) }}

??? toggle "Augur Traits"

{{ include_md("common_text/augur_traits_task.md", indent=8) }}

??? toggle "Extract Clade-defining Mutations"

{{ include_md("common_text/extract_clade_mutations_task.md", indent=8) }}

??? toggle "Augur Clades"

{{ include_md("common_text/augur_clades_task.md", indent=8) }}

??? toggle "Augur Export"

{{ include_md("common_text/augur_export_task.md", indent=8) }}


#### Augur Outputs

!!! dna "Diversity dependent"
    Note that the node & branch coloring by clade or lineage assignment might be dependent on the diversity of your input dataset. This is because the clade assignment is done using the ancestrally reconstructed amino acid or nucleotide changes at the tree nodes rather than a direct sequence-to-reference mutation comparison. You may notice this happening when you get clade/lineage assignments from NextClade when running TheiaCoV workflows, but no clade/lineage assignment on the Augur Auspice tree.

    To get around this issue, you can upload the Augur output file `merged-metadata.tsv` to Auspice that includes the correct clade/lineage assignments to allow for coloring by Clade.

!!! dna "Flu clade assignments"
    Note that for flu, the clade assignment is usually mostly done for the more recent seasonal influenza viruses. Older strains may get an "unassigned" designation for clades. Therefore, it is important to counter check with the NextClade results from TheiaCoV if the lack of clade assignment is due to analyzing older sequences or sequence related.

The `auspice_input_json` is intended to be uploaded to [Auspice](https://auspice.us/) to view the phylogenetic tree. This provides a visualization of the genetic relationships between your set of samples. The `metadata_merged` output can also be uploaded to add context to the phylogenetic visualization. The `combined_assemblies` output can be uploaded to [UShER](https://genome.ucsc.edu/cgi-bin/hgPhyloPlace) to view the samples on a global tree of representative sequences from the public repositories.

The Nextstrain team hosts documentation surrounding the Augur workflow → Auspice visualization here, which details the various components of the Auspice interface: [How data is exported by Augur for visualisation in Auspice](https://docs.nextstrain.org/en/latest/learn/augur-to-auspice.html).

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "Augur"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///

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
