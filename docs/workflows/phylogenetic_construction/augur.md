# Augur

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Name": "[**Augur**](../workflows/phylogenetic_construction/augur.md)"}, columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level", "Dockstore"]) }}

## Augur Workflows

!!! caption "Augur Workflow Diagram"
    ![Augur Workflow Diagram](../../assets/figures/Augur_MAIN.png)

!!! dna inline end "**Helpful resources for epidemiological interpretation**"

    [CDC's COVID-19 Epidemiology Toolkit](https://www.cdc.gov/advanced-molecular-detection/php/training/) is a useful resource for learning more about genomic epidemiology. Some particularly relevant modules to the Augur workflows include:

    - [Module 1.3: How to read a phylogenetic tree](https://www.cdc.gov/advanced-molecular-detection/php/training/module-1-3.html)
    - [Module 3.1: Getting started with Nextstrain](https://www.cdc.gov/advanced-molecular-detection/php/training/module-3-1.html) (which includes Auspice)
    - [Module 3.4: Walking through Nextstrain trees](https://www.cdc.gov/advanced-molecular-detection/php/training/module-3-4.html)

Genomic epidemiology is an important approach to understand and mitigate disease transmission. A critical step in viral genomic epidemiology is generating phylogenetic trees to explore the genetic relationship between viruses on a local, regional, national, or global scale. The Augur workflows enable **viral phylogenetic analysis** by generating phylogenetic trees from genome assemblies and incorporating metadata into a intuitive visual platform via [Auspice](https://auspice.us).

Two workflows are offered: **Augur_Prep_PHB** and **Augur_PHB**. **Augur_Prep_PHB** prepares sample metadata to be visualized alongside the phylogenetic tree produced by **Augur_PHB**. **Augur_PHB** can be run without metadata to produce only a phylogenetic tree. If you want metadata incorporated in your final tree, these workflows must be run sequentially. The final outputs from **Augur_PHB** can be visualized in [Auspice](https://auspice.us), which is the recommended platform. Alternative tree visualization platforms can also be used, though these may not support all metadata features.

### Augur_Prep_PHB

The **Augur_Prep_PHB** workflow was written to prepare individual sample assemblies and their metadata for inclusion in Augur_PHB analysis. The optional metadata inputs include collection date information (in `YYYY-MM-DD` format), clade information (like `nextclade clade` and/or `pango lineage`), and geographical information.

This workflow runs on the sample level, and takes assembly FASTA files and associated metadata formatted in a data table. FASTA files may be generated with one of the TheiaCoV/TheiaViral Characterization workflows and should adhere to quality control guidelines, (e.g. [QC guidelines produced by PHA4GE](https://github.com/pha4ge/pipeline-resources/blob/main/docs/qc-solutions.md)).

!!! dna "**How to prepare metadata**"

    If you are running this workflow on Terra, we recommend carefully preparing metadata in a TSV file that will be uploaded to the same Terra datatable that contains the sample genetic information. An example of a correctly formatted TSV file can be found in [this example](../../assets/files/example-augur-prep-metadata.tsv). A few important considerations are:

    - Please always include the **date information** in `YYYY-MM-DD` format. Other date formats are incompatible with Augur. You can specify unknown dates or month by replacing the respective values by `XX` (e.g. `2013-01-XX` or `2011-XX-XX`), while completely unknown dates can be shown with `20XX-XX-XX` (which does not restrict the sequence to being in the 21st century - they could be earlier). Alternatively, reduced precision format can also be used (e.g. `2018`, `2018-03`).
        - Because Excel will automatically change the date formatting, we recommend not opening or preparing your meta data file in Excel. If the metadata is already in Excel, or you decide to prepare it in Excel, we recommend using another program to correct the dates afterwards (and be caution if you open it in Excel again!).
    - Different levels of **geographical information** can be passed to Augur. A latitude and longitude file is provided by default by Theiagen, which mirrors what can be found [here](https://github.com/nextstrain/augur/blob/master/augur/data/lat_longs.tsv). Just ensure that your spelling matches what is in the file exactly or alternatively provide your own. Augur_Prep supports the following levels:
        - `region` - Lowest-level resolution, used often for continents (e.g.: `europe`, `asia`, `north america`)
        - `country` - Denotes the country where the sample originated (e.g.: `Argentina`, `Japan`, `USA`)
        - `divisions` - Denotes divisions, or states, or sometimes cities, within the country (e.g.: `California`, `Colorado`, `Cork`)
        - `location` - Highest-level resolution, often used for custom latitude and longitude for futher detail on divisions, like cities within states. Just ensure that this level is provided in either the default latitude and longitude file or in a custom one. 
    - Optional **clade** information, such the one assigned by *Nextclade*.
    - Optional **Pangolin lineage** information for SARS-CoV-2 samples.

#### Augur_Prep Inputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "Augur_Prep"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

#### Augur_Prep Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "Augur_Prep"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///

### Augur_PHB

[Augur](https://docs.nextstrain.org/projects/augur/en/stable/) is a bioinformatics toolkit to track evolution from sequence data, ingesting sequences and metadata such as dates and sampling locations, filtering the data, aligning the sequences, infering a tree, and export the results in a format that can be visualized by [Auspice](https://auspice.us/). This is the tool behind [Nextrain's builds](https://nextstrain.org/) available for a large collection of viral organisms.

!!! info "Before getting started"
    **Phylogenetic inference requires careful planning, quality control of sequences, and metadata curation**. You may have to generate phylogenies multiple times by running the Augur_PHB workflow, with several iterations of assessing results and amending inputs, in order to generate a final tree with sufficient diversity and high-quality data of interest. Theiagen's [Introduction to Phylogenetics](../../guides/phylogenetics.md) is one resource to give you the necessary information on the considerations you'll need to have before performing this type of analysis.

**Augur_PHB** takes as input a **set of assembly/consensus files** (FASTA format), an **optional viral organism designation**, and an **optional sample metadata files** (TSV format) that have been formatted via the Augur_Prep_PHB workflow. Augur_PHB runs [Augur](https://docs.nextstrain.org/projects/augur/en/stable/) to generate a phylogenetic tree following the construction of a SNP distance matrix and alignment. Provided metadata will be used to refine the final tree and incorporated into the [Auspice](https://auspice.us/)-formatted tree visual.

#### Augur Inputs

!!! warning "Sample diversity and tree building"
    Before attempting a phylogenetic tree, you must ensure that the input FASTAs meet quality-control metrics. Sets of FASTAs with highly discordant quality metrics may result in the inaccurate inference of genetic relatedness.

    There must be some sequence diversity among the set of input assemblies. If insufficient diversity is present, it may be necessary to add a more divergent sequence to the set of samples to be analyzed.

Some inputs will automatically bypass or trigger modules, such as populating `alignment_fasta`, which bypasses alignment. Clade-defining mutations can be automatically extracted if the "clade_membership" metadata field is provided and the `extract_clade_mutations` optional input is set to `true`.

Any metadata present in the final JSON file for Auspice visualization is determined by what metadata was provided by the user. If the `sample_metadata_tsvs` optional input parameter is **not** provided, the final tree visual will only include the distance tree. If metadata was provided, different metadata fields will trigger different steps: date information will trigger the refinement of the distance tree into a time tree; clade information will be assigned to the tree nodes; geographical information will be represented in the Auspice visual within a map. The following figure illustrates this logic.

!!! caption "Augur Metadata Conditionals"
    ![Augur Metadata Conditionals](../../assets/figures/Augur_Metadata_Conditionals.png){data-description="The metadata and type of tree in the output JSON for Auspice will depend on the metadata that is present in input metadata file. If no metadata file is provided, the output JSON will only contain a distance tree. If date information is present, the distance tree will be replaced by a tree refined by time (time tree). If clade and/or pango lineage (for SARS-CoV-2) information is provided, the tree will display an option to color by lineage. If geographical information is present, a map will load in Auspice using the information provided."}

##### **A Note on Optional Inputs**

!!! warning "Defaults change based on the specified organism input"
    Default values that mimic the Nextstrain builds for the following organisms have been preselected:

    - Flu (`"flu"`), which also requires the following two inputs:
        - `flu_segment` (`"HA"` or `"NA"`)
        - `flu_subtype` (`"H1N1"`, `"H3N2"`, `"Victoria"`, `"Yamagata"`, or `"H5N1"`)
    - RSV-A (`"rsv-a"`)
    - RSV-B (`"rsv-b"`)
    - Mpox (`"mpox"`)
    - SARS-CoV-2 (default; `"sars-cov-2"`) 
 
    View these default parameters in the relevant toggle block below.

{{ include_md("common_text/organism_parameters_wf.md", condition="virus", indent=4) }}

???+ info "Running Augur_PHB on custom organsisms"
    For non-default organisms (listed above), several optional inputs are required to guarantee workflow functionality

    | **Task** | **Input** | **Description** |
    | --- | --- | --- |
    | augur | _reference_fasta_ | Reference sequence in FASTA format. |
    | augur | _reference_genbank_ | Reference sequence in GenBank format. |
    | augur | _organism_ | Name of expected organism. |
    | augur | _min_num_unambig_ | Minimum number of called bases in genome to pass prefilter. | 
    | augur | _lat_longs_tsv_ | Tab-delimited file of geographic location names with corresponding latitude and longitude values. Only necessary if geographical information is in the metadata; must follow [this](https://github.com/nextstrain/augur/blob/master/augur/data/lat_longs.tsv) format |
    | augur | _clades_tsv_ | TSV file containing clade mutation positions in four columns. Only necessary if clade information is in the metadata. |

    In the inputs table below, these fields have both the "Required" and "Optional" tags.

There are **many** optional user inputs. For more information regarding these optional inputs, please view [Nextstrain's detailed documentation on Augur](https://docs.nextstrain.org/projects/augur/en/stable/usage/usage.html)

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "Augur"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Workflow Tasks

For the Augur subcommands, please view the [Nextstrain Augur documentation](https://docs.nextstrain.org/projects/augur/en/stable/usage/usage.html) for more details and explanations.

{{ include_md("common_text/versioning_task.md") }}

{{ include_md("common_text/augur_align_task.md") }}

{{ include_md("common_text/snp_dists_task.md", condition="augur") }}

{{ include_md("common_text/augur_tree_task.md") }}

{{ include_md("common_text/reorder_matrix_task.md") }}

{{ include_md("common_text/augur_refine_task.md") }}

{{ include_md("common_text/augur_ancestral_task.md") }}

{{ include_md("common_text/augur_translate_task.md") }}

{{ include_md("common_text/augur_mutation_context_task.md") }}

{{ include_md("common_text/augur_traits_task.md") }}

{{ include_md("common_text/extract_clade_mutations_task.md") }}

{{ include_md("common_text/augur_clades_task.md") }}

{{ include_md("common_text/augur_export_task.md") }}

#### Augur Outputs

The `auspice_input_json` is intended to be uploaded to [Auspice](https://auspice.us/) to view the resulting phylogenetic tree with the provided metadata. Alternatively, a phylogenetic tree in Newick fromat is also available for visualization in other platforms. The `metadata_merged` output can also be uploaded to either Auspice or a different visualization platform to add further context to the phylogenetic visualization. The `combined_assemblies` output can be uploaded to [UShER](https://genome.ucsc.edu/cgi-bin/hgPhyloPlace) to view the samples on a global tree of representative sequences from public repositories.

The Nextstrain team hosts documentation surrounding the Augur workflow to Auspice visualization here, which details the various components of the Auspice interface: [How data is exported by Augur for visualisation in Auspice](https://docs.nextstrain.org/en/latest/learn/augur-to-auspice.html).

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "Augur"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///

### References

When publishing work using the Augur_PHB workflow, please reference the following:

> Nextstrain: Hadfield J, Megill C, Bell SM, Huddleston J, Potter B, Callender C, Sagulenko P, Bedford T, Neher RA. Nextstrain: real-time tracking of pathogen evolution. Bioinformatics. 2018 Dec 1;34(23):4121-3.
