# Mercury_Prep_N_Batch

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Public Data Sharing](../../workflows_overview/workflows_type.md/#public-data-sharing) | [Viral](../../workflows_overview/workflows_kingdom.md/#viral) | PHB v2.3.0 | Yes | Set-level |

## Mercury_Prep_N_Batch_PHB

Mercury prepares and formats metadata and sequencing files **located in Google Cloud Platform (GCP) buckets** for submission to national & international databases, currently NCBI & GISAID. Mercury was initially developed to ingest read, assembly, and metadata files associated with SARS-CoV-2 amplicon reads from clinical samples and format that data for submission per the [Public Health Alliance for Genomic Epidemiology (PH4GE)'s SARS-CoV-2 Contextual Data Specifications](https://github.com/pha4ge/SARS-CoV-2-Contextual-Data-Specification).

Currently, Mercury supports submission preparation for SARS-CoV-2, mpox, and influenza. These organisms have different metadata requirements, and are submitted to different repositories; the following table lists the repositories for each organism & what is supported in Mercury:

|  | BankIt (NCBI) | BioSample (NCBI) | GenBank (NCBI) | GISAID | SRA (NCBI) |
| --- | --- | --- | --- | --- | --- |
| **`"flu"`** |  | ✓ |  |  | ✓ |
| **`"mpox"`** | ✓ | ✓ |  | ✓ | ✓ |
| **`"sars-cov-2"`** |  | ✓ | ✓ | ✓ | ✓ |

!!! dna "Mercury expects data tables made with TheiaCoV"
    Mercury was designed to work with metadata tables that were partially created after running the TheiaCoV workflows. If you are using a different pipeline, please ensure that the metadata table is formatted correctly. See [this file](https://github.com/theiagen/mercury/blob/main/mercury/Metadata.py) for the hard-coded list of all of the different metadata fields expected for each organism.

### Metadata Formatters

To help users collect all required metadata, we have created the following Excel spreadsheets that can help you collect the necessary metadata and allow for easy upload of this metadata into your Terra data tables:

??? toggle "**_For flu_**"

    [Flu Metadata Formatter](../../assets/metadata_formatters/Terra_2_NCBI-PATHOGEN-metadata-2024-04-30.xlsx)

    Flu uses the same metadata formatter as the Terra_2_NCBI Pathogen BioSample package.

    If neither `strain` nor `isolate` are found in the Terra data table, Mercury will automatically generate an isolate, using the following format 
    `ABRicate flu type / State / sample name / year (ABRicate flu subtype)`. Example: `A/California/Sample-01/2024 (H1N1)`

    The ABRicate flu type and subtype (`abricate_flu_type` and `abricate_flu_subtype` columns) are extracted from your table, and are required to generate the isolate field if it is not provided.

??? toggle "**_For mpox_**"

    [Mpox Metadata Formatter](../../assets/metadata_formatters/Mercury_Prep_N_Batch_MPXV_Metadata_Formatter_2022_12_23.xlsx)

??? toggle "**_For sars-cov-2_**"

    [SARS-CoV-2 Metadata Formatter](../../assets/metadata_formatters/Mercury_Prep_N_Batch_SC2_Metadata_Formatter_2023_05_22.xlsx)

!!! dna "Usage on Terra"

    ##### Usage on Terra {% raw %} {#usage-on-terra} {% endraw %}

    **A note on the `using_clearlabs_data` & `using_reads_dehosted` optional input parameters**

    The `using_clearlabs_data` and `using_reads_dehosted` arguments change the default values for the `read1_column_name`, `assembly_fasta_column_name`, and `assembly_mean_coverage_column_name` metadata columns. The default values are shown in the table below in addition to what they are changed to depending on what arguments are used.

    | Variable | Default Value | with `using_clearlabs_data` | with `using_reads_dehosted` | with both  `using_clearlabs_data` **_and_** `using_reads_dehosted` |
    | --- | --- | --- | --- | --- |
    | `read1_column_name` | `"read1_dehosted"` | `"clearlabs_fastq_gz"` | `"reads_dehosted"` | `"reads_dehosted"` |
    | `assembly_fasta_column_name` | `"assembly_fasta"` | `"clearlabs_fasta"` | `"assembly_fasta"` | `"clearlabs_fasta"` |
    | `assembly_mean_coverage_column_name` | `"assembly_mean_coverage"` | `"clearlabs_sequencing_depth"` | `"assembly_mean_coverage"` | `"clearlabs_sequencing_depth"` |

### Inputs

!!! tip "Use the sample table for the `terra_table_name` input"
    Make sure your entry for `terra_table_name` is for the _sample_ table! While the root entity needs to be the set table, the input value for `terra_table_name` should be the sample table.

This workflow runs on the set-level.

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "Mercury_Prep_N_Batch"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "Mercury_Prep_N_Batch"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///

???+ toggle "An example excluded_samples.tsv file"

    ##### An example excluded_samples.tsv file {% raw %} {#example-excluded-samples} {% endraw %}

    Due to the nature of tsv files, it may be easier to download and open this file in Excel. 
    
    [example_excluded_samples.tsv](../../assets/files/example_excluded_samples.tsv)
    
    ```
    Samples excluded for quality thresholds:
    sample_name message 
    sample2 VADR skipped due to poor assembly
    sample3 VADR number alerts too high: 3 greater than limit of 0
    sample4 Number of Ns was too high: 10000 greater than limit of 5000
    
    Samples excluded for missing required metadata (will have empty values in indicated columns):
    tablename_id    organism    country library_layout
    sample5         paired
    sample6 SARS-CoV-2  USA
    ```
    
    This example informs the user that samples 2-4 were excluded for quality reasons (the exact reason is listed in the `message` column), and that samples 5 and 6 were excluded because they were missing required metadata fields (sample5 was missing the `organism` and `country` fields, and sample6 was missing the `library_layout` field).

## Usage outside of Terra

This tool can also be used on the command-line. Please see [the Mercury GitHub](https://github.com/theiagen/mercury) for more information on how to run Mercury with a Docker image or in your local command-line environment.
