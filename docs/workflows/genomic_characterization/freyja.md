# Freyja Workflow Series

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Genomic Characterization](../../workflows_overview/workflows_type.md/#genomic-characterization) | [Viral](../../workflows_overview/workflows_kingdom.md/#viral) | PHB v3.0.0 | Yes | Sample-level, Set-level |

## Freyja Overview

[Freyja](https://github.com/andersen-lab/Freyja) is a tool for analysing viral mixed sample genomic sequencing data. Developed by Joshua Levy from the [Andersen Lab](https://andersen-lab.com/), it performs two main steps:

1. **Variant Frequency Estimation:** Freyja calculates the frequencies of single nucleotide variants (SNVs) in the genomic sequencing data.
2. **Depth-Weighted Demixing:** It separates mixed populations of viral subtypes using a depth-weighted statistical approach, estimating the proportional abundance of each subtype in the sample based on the frequencies of subtype-defining variants.

Additional post-processing steps can produce visualizations of aggregated samples.

!!! dna "Wastewater and more"
    The typical use case of Freyja is to **analyze mixed SARS-CoV-2 samples** from a sequencing dataset, most often **wastewater**, but the tool is not limited to this context. With the appropriate reference genomes and barcode files, Freyja can be adapted for other pathogens, including MPXV, Influenza, RSV, and Measles.

    !!! warning "Default Values"
        The defaults included in the Freyja workflows reflect this use case but **can be adjusted for other pathogens**. See the [**Running Freyja on other pathogens**](freyja.md#running-freyja-on-other-pathogens) section for more information. Please be aware this is an ==_experimental feature_== and we cannot guarantee complete functionality at this time.

!!! caption "Figure 1: Workflow diagram for Freyja Suite of workflows"
    ##### Figure 1 { #figure1 }
    ![**Figure 1: Workflow diagram for Freyja Suite of workflows.**](../../assets/figures/Freyja_Suite.png){width=100%}

    Depending on the type of data (Illumina or Oxford Nanopore), the Read QC and Filtering steps, as well as the Read Alignment steps use different software. The user can specify if the barcodes and lineages file should be updated with `freyja update` before running Freyja or if bootstrapping is to be performed with `freyja boot`.

Four workflows have been created that perform different parts of Freyja:

- [**Freyja_Update_PHB**](freyja.md#freyja_update)
- [**Freyja_FASTQ_PHB**](freyja.md#freyja_fastq)
- [**Freyja_Plot_PHB**](freyja.md#freyja_plot)
- [**Freyja_Dashboard_PHB**](freyja.md#freyja_dashboard)

The main workflow is [**Freyja_FASTQ_PHB**](freyja.md#freyja_fastq) ([Figure 1](freyja.md#figure1)). Depending on the type of input data (Illumina paired-end, Illumina single-end or ONT), it runs various QC modules before aligning the sample with either [BWA](https://github.com/lh3/bwa) (Illumina) or [minimap2](https://github.com/lh3/minimap2) (ONT) to the provided reference file, followed by iVar for primer trimming. After the preprocessing is completed, [Freyja](https://github.com/andersen-lab/Freyja) is run to generate relative lineage abundances (demix) from the sample. Optional bootstrapping may be performed.

!!! dna "Data Compatability"

    The **Freyja_FASTQ_PHB workflow** is compatible with the following input data types:

        - Ilumina Single-End
        - Illumina Paired-End
        - Oxford Nanopore

[**Freyja_Update_PHB**](freyja.md#freyja_update) will copy the **SARS-CoV-2** reference files that can then be used as input for the [Freyja_FASTQ_PHB](freyja.md#freyja_fastq) workflow.

Two options are available to visualize the Freyja results: [**Freyja_Plot_PHB**](freyja.md#freyja_plot) and [**Freyja_Dashboard_PHB**](freyja.md#freyja_dashboard). [Freyja_Plot_PHB](freyja.md#freyja_plot) aggregates multiple samples using output from [Freyja_FASTQ_PHB](freyja.md#freyja_fastq) to generate a plot that shows fractional abundance estimates for all samples. including the option to plot sample collection date information. Alternatively, [**Freyja_Dashboard_PHB**](freyja.md#freyja_dashboard) aggregates multiple samples using output from [Freyja_FASTQ_PHB](freyja.md#freyja_fastq) to generate an interactive visualization. This workflow requires an additional input field called viral load, which is the number of viral copies per liter.

### Freyja, Sequencing Platforms and Data Quality

The choice of sequencing platform and the quality of the data directly influence Freyja's performance. High-accuracy platforms like Illumina provide reliable SNV detection, enhancing the precision of lineage abundance estimates. In contrast, platforms with higher error rates, such as Nanopore, whilst it has improved greatly in the recent years, may introduce uncertainties in variant calling, affecting the deconvolution process. Sequencing depth requirements will increase as the quality of the sequencing data decreases. A rational target depth is 100X coverage for sequencing data with Q-scores in the range of 25-30.

Additionally, inadequate sequencing depth can hinder Freyja's ability to differentiate between lineages, leading to potential misestimations. Sequencing depth requirements will increase with the complexity of the sample composition and the diversity of lineages present. For samples containing multiple closely related lineages, higher sequencing depth is necessary to resolve subtle differences in genetic variation and accurately estimate lineage abundances. This is particularly important for pathogens with high mutation rates or a large number of cocirculating lineages, such as influenza, where distinguishing between lineages relies on detecting specific single nucleotide variants (SNVs) with high confidence.

## Freyja Workflows

### Freyja_Update_PHB {% raw %} {#freyja_update} {% endraw %}

This workflow will copy the **SARS-CoV-2 reference files** (`curated_lineages.json` and `usher_barcodes.feather`) from [the source repository](https://github.com/andersen-lab/Freyja/tree/main/freyja/data) to a user-specific Google Cloud Storage (GCP) location (often a [Terra.bio](http://Terra.bio) workspace-associated bucket). These files can then be used as input for the [Freyja_FASTQ_PHB workflow](freyja.md#freyja_fastq).

!!! warning "Warning"

    This workflow is compatible only with **SARS-CoV-2 reference files**! To download reference files for other organisms please see the following repository: [Freyja Barcodes](https://github.com/gp201/Freyja-barcodes).

    More information is available in the [**Running Freyja on other pathogens**](freyja.md#running-freyja-on-other-pathogens) section.

#### Inputs

We recommend running this workflow with **"Run inputs defined by file paths"** selected since no information from a Terra data table is actually being used. We also recommend turning off call caching so new information is retrieved every time.

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/input_tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="Freyja_Update", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

#### Outputs

This workflow does not produce any outputs that appear in a Terra data table. The reference files will appear at the location specified with the `gcp_uri` input variable.

### Freyja_FASTQ_PHB {% raw %} {#freyja_fastq} {% endraw %}

Freyja measures SNV frequency and sequencing depth at each position in the genome to return an estimate of the true lineage abundances in the sample. The method uses lineage-defining "barcodes" that, for SARS-CoV-2, are derived from the UShER global phylogenetic tree as a base set for demixing. **Freyja_FASTQ_PHB** returns as output a TSV file that includes the lineages present and their corresponding abundances, along with other values.

The Freyja_FASTQ_PHB workflow is compatible with the multiple input data types: Ilumina Single-End, Illumina Paired-End and Oxford Nanopore. Depending on the type of input data, different input values are used.

**Table 1:** Freyja_FASTQ_PHB input configuration for different types of input data.

| Table Columns | Illumina Paired-End | Illumina Single-End | Oxford Nanopore |
| --- | --- | --- | --- |
| **read1** | ✅ | ✅ | ✅ |
| **read2** | ✅ | ❌ | ❌ |
| **ont** | `false` | `false` | `true` |

#### Inputs

This workflow runs on the sample level.

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/input_tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="Freyja_FASTQ", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

#### Analysis Tasks

??? task "`read_QC_trim_pe` Details"

    `read_QC_trim_pe` {% raw %} {#read_QC_trim_pe} {% endraw %}

    This task runs a sub-workflow that gathers basic QC information, trimming (either with trimmomatic or fastp), human read scrubbing, and taxonomic identification (Kraken2). Optional parameters do not need to be modified. For information regarding the individual tasks performed during this, please visit the [TheiaCoV documentation](../genomic_characterization/theiacov.md).

    !!! techdetails "Read_QC_Trim_PE Technical Details"

        |  | Links |
        | --- | --- |
        | Task | [wf_read_QC_trim_pe.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_read_QC_trim_pe.wdl) |

??? task "`read_QC_trim_se` Details"

    `read_QC_trim_se` {% raw %} {#read_QC_trim_se} {% endraw %}

    This task runs a sub-workflow that gathers basic QC information, trimming (either with trimmomatic or fastp), human read scrubbing, and taxonomic identification (Kraken2). Optional parameters do not need to be modified. For information regarding the individual tasks performed during this, please visit the [TheiaCoV documentation](../genomic_characterization/theiacov.md).

    !!! techdetails "Read_QC_Trim_SE Technical Details"

        |  | Links |
        | --- | --- |
        | Task | [wf_read_QC_trim_se.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_read_QC_trim_se.wdl) |

??? task "`read_QC_trim_ont` Details"

    `read_QC_trim_ont` {% raw %} {#read_QC_trim_ont} {% endraw %}

    This task runs a sub-workflow that gathers basic QC information, trimming (nanoplot), human read scrubbing, and taxonomic identification (Kraken2). Optional parameters do not need to be modified. For information regarding the individual tasks performed during this, please visit the [TheiaCoV documentation](../genomic_characterization/theiacov.md).

    !!! techdetails "Read_QC_Trim_ONT Technical Details"
    
        |  | Links |
        | --- | --- |
        | Task | [wf_read_QC_trim_ont.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_read_QC_trim_ont.wdl) |

??? task "`bwa` Details"

    `bwa` {% raw %} {#bwa} {% endraw %}

    This task aligns the cleaned short reads (Illumina) to the reference genome provided by the user.

    !!! techdetails "BWA Technical Details"
    
        |  | Links |
        | --- | --- |
        | Task | [task_bwa.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/alignment/task_bwa.wdl) |
        | Software Source Code | https://github.com/lh3/bwa |
        | Software Documentation | https://bio-bwa.sourceforge.net/ |
        | Original Publication(s) | [Fast and accurate short read alignment with Burrows-Wheeler transform](https://doi.org/10.1093/bioinformatics/btp324) |

??? task "`minimap2` Details"

    `minimap2` {% raw %} {#minimap2} {% endraw %}

    This task aligns the cleaned long reads (Oxford Nanopore) to the reference genome provided by the user.

    !!! techdetails "Minimap2 Technical Details"
    
        |  | Links |
        | --- | --- |
        | Task | [task_minimap2.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/alignment/task_minimap2.wdl) |
        | Software Source Code | https://github.com/lh3/minimap2 |
        | Software Documentation | https://lh3.github.io/minimap2/ |
        | Original Publication(s) | [Minimap2: pairwise alignment for nucleotide sequences](https://doi.org/10.1093/bioinformatics/bty191) |

??? task "`primer_trim` Details"

    `primer_trim` {% raw %} {#primer_trim} {% endraw %}

    This task trims the primer sequences from the aligned bam file with iVar. The optional input, `keep_noprimer_reads`, does not have to be modified.

    !!! techdetails "Primer Trim Technical Details"
    
        |  | Links |
        | --- | --- |
        | Task | [task_ivar_primer_trim.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_ivar_primer_trim.wdl) |
        | Software Source Code | https://github.com/andersen-lab/ivar |
        | Software Documentation | https://andersen-lab.github.io/ivar/html/manualpage.html |
        | Original Publication(s) | [An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar](https://doi.org/10.1186/s13059-018-1618-7) |

??? task "`freyja` Details"

    `freyja` {% raw %} {#freyja} {% endraw %}

    The Freyja task will call variants and capture sequencing depth information to identify the relative abundance of lineages present. Optionally, if `bootstrap` is set to true, bootstrapping will be performed. After the optional bootstrapping step, the variants are demixed.

    !!! techdetails "Freyja Technical Details"

        |  | Links |
        | --- | --- |
        | Task | [task_freyja_one_sample.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/freyja/task_freyja.wdl) |
        | Software Source Code | https://github.com/andersen-lab/Freyja |
        | Software Documentation | https://andersen-lab.github.io/Freyja/index.html# |

#### Outputs

The main output file used in subsequent Freyja workflows is found under the `freyja_demixed` column. This TSV file takes on the following format:

|  | sample name |
| --- | --- |
| summarized | [('Delta', 0.65), ('Other', 0.25), ('Alpha', 0.1')] |
| lineages | ['B.1.617.2' 'B.1.2' 'AY.6' 'Q.3'] |
| abundances | "[0.5 0.25 0.15 0.1]" |
| resid | 3.14159 |
| coverage | 95.8 |

- The `summarized` array denotes a sum of all lineage abundances in a particular WHO designation (i.e. B.1.617.2 and AY.6 abundances are summed in the above example), otherwise they are grouped into "Other".
- The `lineage` array lists the identified lineages in descending order
- The `abundances` array contains the corresponding abundances estimates.
- The value of `resid` corresponds to the residual of the weighted least absolute deviation problem used to estimate lineage abundances.
- The `coverage` value provides the 10x coverage estimate (percent of sites with 10 or greater reads)

!!! tip "Click "Ignore empty outputs""
    When running the Freyja_FASTQ_PHB workflow, it is recommended to select the "Ignore empty outputs" option in the Terra UI. This will hide the output columns that will not be generated for your input data type.

<div class="searchable-table" markdown="1">

| **Variable** | **Type** | **Description** | **Input Data Type** |
|---|---|---|---|
| aligned_bai | File | Index companion file to the bam file generated during the consensus assembly process | ONT, PE, SE |
| aligned_bam | File | Primer-trimmed BAM file; generated during consensus assembly process | ONT, PE, SE |
| alignment_method | String | The method used to generate the alignment | ONT, PE, SE |
| bbduk_docker | String | Docker image used to run BBDuk | PE, SE |
| bwa_version | String | Version of BWA used to map read data to the reference genome | PE, SE |
| fastp_html_report | File | The HTML report made with fastp | PE, SE |
| fastp_version | String | Version of fastp software used | PE, SE |
| fastq_scan_clean1_json | File | JSON file output from `fastq-scan` containing summary stats about clean forward read quality and length | PE, SE |
| fastq_scan_clean2_json | File | JSON file output from `fastq-scan` containing summary stats about clean reverse read quality and length | PE |
| fastq_scan_num_reads_clean_pairs | String | Number of clean read pairs | PE |
| fastq_scan_num_reads_clean1 | Int | Number of clean forward reads | PE, SE |
| fastq_scan_num_reads_clean2 | Int | Number of clean reverse reads | PE |
| fastq_scan_num_reads_raw_pairs | String | Number of raw read pairs | PE |
| fastq_scan_num_reads_raw1 | Int | Number of raw forward reads | PE, SE |
| fastq_scan_num_reads_raw2 | Int | Number of raw reverse reads | PE |
| fastq_scan_raw1_json | File | JSON file output from `fastq-scan` containing summary stats about raw forward read quality and length | PE, SE |
| fastq_scan_raw2_json | File | JSON file output from `fastq-scan` containing summary stats about raw reverse read quality and length | PE |
| fastq_scan_version | String | Version of fastq_scan used for read QC analysis | PE, SE |
| fastqc_clean1_html | File | Graphical visualization of clean forward read quality from fastqc to open in an internet browser | PE, SE |
| fastqc_clean2_html | File | Graphical visualization of clean reverse read quality from fastqc to open in an internet browser | PE |
| fastqc_docker | String | Docker container used for fastqc | PE, SE |
| fastqc_num_reads_clean_pairs | String | Number of read pairs after cleaning by fastqc | PE |
| fastqc_num_reads_clean1 | Int | Number of forward reads after cleaning by fastqc | PE, SE |
| fastqc_num_reads_clean2 | Int | Number of reverse reads after cleaning by fastqc | PE |
| fastqc_num_reads_raw_pairs | String | Number of input read pairs by fastqc | PE |
| fastqc_num_reads_raw1 | Int | Number of input forward reads by fastqc | PE, SE |
| fastqc_num_reads_raw2 | Int | Number of input reverse reads by fastqc | PE |
| fastqc_raw1_html | File | Graphical visualization of raw forward read quality from fastqc to open in an internet browser | PE, SE |
| fastqc_raw2_html | File | Graphical visualization of raw reverse read qualityfrom fastqc to open in an internet browser | PE |
| fastqc_version | String | Version of fastqc software used | PE, SE |
| freyja_abundances | String | Abundances estimates identified by Freyja and parsed from freyja_demixed file | ONT, PE, SE |
| freyja_sc2_barcode_file | File | Barcode feather file used for SARS-CoV-2. Can be the one provided as input or downloaded by Freyja if update_db is true | ONT, PE, SE |
| freyja_barcode_version | String | Name of barcode file used, or the date if update_db is true | ONT, PE, SE |
| freyja_bootstrap_lineages | File | A CSV that contains the 0.025, 0.05, 0.25, 0.5 (median), 0.75, 0.95, and 0.975 percentiles for each lineage | ONT, PE, SE |
| freyja_bootstrap_lineages_pdf | File | A boxplot of the bootstrap lineages CSV file | ONT, PE, SE |
| freyja_bootstrap_summary | File | A CSV that contains the 0.025, 0.05, 0.25, 0.5 (median), 0.75, 0.95, and 0.975 percentiles for each WHO designated VOI/VOC | ONT, PE, SE |
| freyja_bootstrap_summary_pdf | File | A boxplot of the bootstrap summary CSV file | ONT, PE, SE |
| freyja_coverage | Float | Coverage identified by Freyja and parsed from freyja_demixed file | ONT, PE, SE |
| freyja_demixed | File | The main output TSV; see the section directly above this table for an explanation | ONT, PE, SE |
| freyja_demixed_parsed | File | Parsed freyja_demixed file, containing the same information, for easy result concatenation | ONT, PE, SE |
| freyja_depths | File | A TSV listing the depth of every position | ONT, PE, SE |
| freyja_fastq_wf_analysis_date | String | Date of analysis | ONT, PE, SE |
| freyja_fastq_wf_version | String | The version of the Public Health Bioinformatics (PHB) repository used | ONT, PE, SE |
| freyja_sc2_lineage_metadata_file | File | Lineage metadata JSON file used for SARS-CoV-2. Can be the one provided as input or downloaded by Freyja if update_db is true | ONT, PE, SE |
| freyja_lineages | String |  Lineages in descending order identified by Freyja and parsed from freyja_demixed file | ONT, PE, SE |
| freyja_metadata_version | String | Name of lineage metadata file used, or the date if update_db is true | ONT, PE, SE |
| freyja_resid | String | Residual of the weighted least absolute deviation problem used to estimate lineage abundances identified by Freyja and parsed from freyja_demixed file | ONT, PE, SE |
| freyja_summarized | String | Sum of all lineage abundances in a particular WHO designation identified by Freyja and parsed from freyja_demixed file | ONT, PE, SE |
| freyja_variants | File | The TSV file containing the variants identified by Freyja | ONT, PE, SE |
| freyja_version | String | version of Freyja used | ONT, PE, SE |
| ivar_version_primtrim | String | Version of iVar for running the iVar trim command | ONT, PE, SE |
| kraken_human | Float | Percent of human read data detected using the Kraken2 software | ONT, PE, SE |
| kraken_human_dehosted | Float | Percent of human read data detected using the Kraken2 software after host removal | ONT, PE, SE |
| kraken_report | File | Full Kraken report | ONT, PE, SE |
| kraken_report_dehosted | File | Full Kraken report after host removal | ONT, PE, SE |
| kraken_sc2 | String | Percent of SARS-CoV-2 read data detected using the Kraken2 software | ONT, PE, SE |
| kraken_sc2_dehosted | String | Percent of SARS-CoV-2 read data detected using the Kraken2 software after host removal | ONT, PE, SE |
| kraken_version | String | Version of Kraken software used | ONT, PE, SE |
| minimap2_docker | String | Docker image used to run minimap2 | ONT |
| minimap2_version | String | Version of minimap2 used | ONT |
| nanoplot_html_clean | File | Clean read file | ONT |
| nanoplot_html_raw | File | Raw read file | ONT |
| nanoplot_num_reads_clean1 | Int | Number of clean reads for the forward-facing file | ONT |
| nanoplot_num_reads_raw1 | Int | Number of reads for the forward-facing file | ONT |
| nanoplot_r1_est_coverage_clean | Float | Estimated coverage on the clean reads by nanoplot | ONT |
| nanoplot_r1_est_coverage_raw | Float | Estimated coverage on the raw reads by nanoplot | ONT |
| nanoplot_r1_mean_q_clean | Float | Mean quality score of clean forward reads | ONT |
| nanoplot_r1_mean_q_raw | Float | Mean quality score of raw forward reads | ONT |
| nanoplot_r1_mean_readlength_clean | Float | Mean read length of clean forward reads | ONT |
| nanoplot_r1_mean_readlength_raw | Float | Mean read length of raw forward reads | ONT |
| nanoplot_r1_median_q_clean | Float | Median quality score of clean forward reads | ONT |
| nanoplot_r1_median_q_raw | Float | Median quality score of raw forward reads | ONT |
| nanoplot_r1_median_readlength_clean | Float | Median read length of clean forward reads | ONT |
| nanoplot_r1_median_readlength_raw | Float | Median read length of raw forward reads | ONT |
| nanoplot_r1_n50_clean | Float | N50 of clean forward reads | ONT |
| nanoplot_r1_n50_raw | Float | N50 of raw forward reads | ONT |
| nanoplot_r1_stdev_readlength_clean | Float | Standard deviation read length of clean forward reads | ONT |
| nanoplot_r1_stdev_readlength_raw | Float | Standard deviation read length of raw forward reads | ONT |
| nanoplot_tsv_clean | File | Output TSV file created by nanoplot | ONT |
| nanoplot_tsv_raw | File | Output TSV file created by nanoplot | ONT |
| nanoq_version | String | Version of nanoq used in analysis | ONT |
| primer_bed_name | String | Name of the primer bed file used for primer trimming | ONT, PE, SE |
| primer_trimmed_read_percent | Float | Percentage of read data with primers trimmed as determined by iVar trim | ONT, PE, SE |
| read1_clean | File | Forward read file after quality trimming and adapter removal | ONT, PE, SE |
| read1_dehosted | File | Dehosted forward reads | ONT, PE, SE |
| read2_clean | File | Reverse read file after quality trimming and adapter removal | PE |
| read2_dehosted | File | Dehosted reverse reads | PE |
| samtools_version | String | The version of SAMtools used to sort and index the alignment file | ONT, PE, SE |
| samtools_version_primtrim | String | The version of SAMtools used to create the pileup before running iVar trim | ONT, PE, SE |
| trimmomatic_docker | String | Docker container for Trimmomatic | PE, SE |
| trimmomatic_version | String | The version of Trimmomatic used | PE, SE |

</div>

### Freyja_Plot_PHB {% raw %} {#freyja_plot} {% endraw %}

This workflow visualizes aggregated freyja_demixed output files produced by [Freyja_FASTQ_PHB](freyja.md#freyja_fastq) in a single plot (pdf format) which provides fractional abundance estimates for all aggregated samples.

Options exist to provide lineage-specific breakdowns and/or sample collection time information.

#### Inputs

This workflow runs on the set level.


/// html | div[class="searchable-table"]

{{ input_table("docs/assets/input_tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="Freyja_Plot", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

#### Analysis Tasks

??? task "`freyja_plot_task` Details"

    `freyja_plot_task` {% raw %} {#freyja_plot_task} {% endraw %}

    This task will aggregate multiple samples together, and then creates a plot. Several optional inputs dictate the plot appearance (see each variable's description for more information).

    !!! techdetails "Freyja Plot Technical Details"

        |  | Links |
        | --- | --- |
        | Task | [wf_freyja_plot.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/freyja/task_freyja_plot.wdl) |
        | Software Source Code | https://github.com/andersen-lab/Freyja |
        | Software Documentation | https://github.com/andersen-lab/Freyja |

#### Outputs

| **Variable** | **Type** | **Description** |
|---|---|---|
| **freyja_demixed_aggregate** | File | A TSV file that summarizes the `freyja_demixed` otuputs for all samples |
| **freyja_plot** | File | A PDF of the plot produced by the workflow |
| **freyja_plot_metadata** | File | The metadata used to create the plot |
| **freyja_plot_version** | String | The version of Freyja used |
| **freyja_plot_wf_analysis_date** | String | The date of analysis |
| **freyja_plot_wf_version** | String | The version of the Public Health Bioinformatics (PHB) repository used |

### Freyja_Dashboard_PHB {% raw %} {#freyja_dashboard} {% endraw %}

This workflow creates a group of interactive visualizations based off of the aggregated freyja_demixed output files produced by [Freyja_FASTQ_PHB](freyja.md#freyja_fastq) called a "dashboard". Creating this dashboard requires knowing the viral load of your samples (viral copies/litre).

!!! warning

    This dashboard is not "live" — that is, you must rerun the workflow every time you want new data to be included in the visualizations.

#### Inputs

This workflow runs on the set level.


/// html | div[class="searchable-table"]

{{ input_table("docs/assets/input_tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="Freyja_Dashboard", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

#### Analysis Tasks

??? task "`freyja_dashboard_task` Details"

    This task will aggregate multiple samples together, and then create an interactive HTML visualization. Several optional inputs dictate the dashboard appearance (see each variable's description for more information).

    !!! techdetails "Freyja Dashboard Technical Details"
    
        |  | Links |
        | --- | --- |
        | Task | [wf_freyja_dashboard.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/freyja/task_freyja_dashboard.wdl) |
        | Software Source Code | https://github.com/andersen-lab/Freyja |
        | Software Documentation | https://github.com/andersen-lab/Freyja |

#### Outputs

| **Variable** | **Type** | **Description** |
|---|---|---|
| freyja_dashboard | File | The HTML file of the dashboard created |
| freyja_dashboard_metadata | File | The metadata used to create the dashboard |
| freyja_dashboard_version | String | The version of Freyja used |
| freyja_dashboard_wf_analysis_date | String | The date of analysis |
| freyja_dashboard_wf_version | String | The version of the Public Health Bioinformatics (PHB) repository used |
| freyja_demixed_aggregate | File | A TSV file that summarizes the `freyja_demixed` outputs for all samples |

## Running Freyja on other pathogens

!!! warning "Experimental Feature"
    Please be aware this is an _experimental feature_ and we cannot guarantee complete functionality at this time.

The main requirement to run Freyja on other pathogens is **the existence of a barcode file for your pathogen of interest**. Currently, barcodes exist for the following organisms:

- SARS-CoV-2 (default)
- FLU-B-VIC
- H1N1
- H3N2
- H5Nx-cattle
- H5NX
- MEASLESN450
- MEASLESgenome
- MPX
- RSVa
- RSVb

!!! dna "Freyja barcodes for other pathogens"

    Data for various pathogens can be found in the following repository: [Freyja Barcodes](https://github.com/gp201/Freyja-barcodes)

    Folders are organized by pathogen, with each subfolder named after the date the barcode was generated, using the format YYYY-MM-DD, as well as a "latest" folder. Barcode files are named `barcode.csv`, and reference genome files are named `reference.fasta`.

There are two ways to run [**Freyja_FASTQ_PHB**](freyja.md#freyja_fastq) for non-SARS-CoV-2 organisms:
- Using the `freyja_pathogen` optional input (limited set of allowable organisms)
- Providing the appropriate barcode file through the `freyja_barcodes` optional input (any organism for which barcodes are supplied)

### Using the `freyja_pathogen` flag

When using the `freyja_pathogen` flag, the user must set the optional `update_db` flag to _true_, so that the latest version of the barcode file is automatically downloaded by Freyja. 

!!! caption "Figure 2:  Optional input for Freyja_FASTQ_PHB to provide the pathogen to be used by Freyja"
    ##### Figure 2 { #figure2 }
    ![**Figure 2:  Optional input for Freyja_FASTQ_PHB to provide the pathogen to be used by Freyja.**](../../assets/figures/Freyja_figure2.png)

Allowed options:

- SARS-CoV-2 (default)
- MPXV
- H1N1pdm
- H5NX
- FLU-B-VIC
- MEASLESN450
- MEASLES
- RSVa
- RSVb

!!! warning

    The `freyja_pathogen` flag is not used if a barcodes file is provided. This means that this option is ignored if a barcode file is provided through `freyja_barcodes`.

### Providing the appropriate barcode file

The appropriate barcode file for your organism of interest and reference sequence need to be downloaded and uploaded to your [Terra.bio](http://Terra.bio) workspace. When running [**Freyja_FASTQ_PHB**](freyja.md#freyja_fastq), the appropriate reference and barcodes file need to be passed as inputs. The first is a required input and will show up at the top of the workflows inputs page on [Terra.bio](http://Terra.bio) ([Figure 3](freyja.md/#figure3)).

!!! caption "Figure 3:  Required input for Freyja_FASTQ_PHB to provide the reference genome to be used by Freyja"
    ##### Figure 3 {% raw %} {#figure3} {% endraw %}
    ![**Figure 3:  Required input for Freyja_FASTQ_PHB to provide the reference genome to be used by Freyja.**](../../assets/figures/Freyja_figure3.png)

The barcodes file can be passed directly to Freyja by the `freyja_barcodes` optional input ([Figure 4](freyja.md/#figure4)).

!!! caption "Figure 4: Optional input for Freyja_FASTQ_PHB to provide the barcodes file to be used by Freyja"
    ##### Figure 4 {% raw %} {#figure4} {% endraw %}
    ![**Figure 4: Optional input for Freyja_FASTQ_PHB to provide the barcodes file to be used by Freyja.**](../../assets/figures/Freyja_figure4.png)

## References

If you use any of the Freyja workflows, please cite:

> Karthikeyan, S., Levy, J.I., De Hoff, P. *et al.* Wastewater sequencing reveals early cryptic SARS-CoV-2 variant transmission. *Nature* **609**, 101–108 (2022). <https://doi.org/10.1038/s41586-022-05049-6>
<!-- -->
> Freyja source code can be found at <https://github.com/andersen-lab/Freyja>
<!-- -->
> Freyja barcodes (non-SARS-CoV-2): <https://github.com/gp201/Freyja-barcodes>
