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

    ##### Usage on Terra {#usage-on-terra}

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

<div class="searchable-table" markdown="1">

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| mercury_prep_n_batch | **gcp_bucket_uri** | String | Google bucket where your SRA reads will be temporarily stored before transferring to SRA. Example: "gs://theiagen_sra_transfer" |  | Required |
| mercury_prep_n_batch | **sample_names** | Array[String] | The samples you want to submit |  | Required |
| mercury_prep_n_batch | **terra_project_name** | String | The name of your Terra project. You can find this information in the URL of the webpage of your Terra dashboard. For example, if your URL contains `#workspaces/example/my_workspace/` then your project name is `example` |  | Required |
| mercury_prep_n_batch | **terra_table_name** | String | The name of the Terra table where your **samples** can be found. Do not include the `entity:` prefix, the `_id` suffix, or the `_set_id` suffix, just the name of the sample-level data table as listed in the sidebar on lefthand side of the Terra Data tab. |  | Required |
| mercury_prep_n_batch | **terra_workspace_name** | String | The name of your Terra workspace where your samples can be found. For example, if your URL contains #workspaces/example/my_workspace/ then your project name is my_workspace |  | Required |
| download_terra_table | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional |
| download_terra_table | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 10 | Optional |
| download_terra_table | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-06-21 | Optional |
| download_terra_table | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional |
| mercury | **amplicon_primer_scheme** | String | Populate to overwrite `amplicon_primer_scheme` column |  | Optional |
| mercury | **amplicon_size** | String | Populate to overwrite `amplicon_size` column |  | Optional |
| mercury | **authors** | String | Populate to overwrite `authors` column |  | Optional |
| mercury | **bioproject_accession** | String | Populate to overwrite `bioproject_accession` column |  | Optional |
| mercury | **continent** | String | Populate to overwrite `continent` column |  | Optional |
| mercury | **country** | String | Populate to overwrite `country` column |  | Optional |
| mercury | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| mercury | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| mercury | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/mercury:1.1.0 | Optional |
| mercury | **gisaid_submitter** | String | Populate to overwrite `gisaid_submitter` column |  | Optional |
| mercury | **host_disease** | String | Populate to overwrite `host_disease` column |  | Optional |
| mercury | **instrument_model** | String | Populate to overwrite `instrument_model` column |  | Optional |
| mercury | **isolation_source** | String | Populate to overwrite `isolation_source` column |  | Optional |
| mercury | **library_layout** | String | Populate to overwrite `library_layout` column |  | Optional |
| mercury | **library_selection** | String | Populate to overwrite `library_selection` column |  | Optional |
| mercury | **library_source** | String | Populate to overwrite `library_source` column |  | Optional |
| mercury | **library_strategy** | String | Populate to overwrite `library_strategy` column |  | Optional |
| mercury | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| mercury | **metadata_organism** | String | Populate to overwrite `organism` column |  | Optional |
| mercury | **number_N_threshold** | Int | Only for "sars-cov-2" submissions; used to filter out any samples that contain more than the indicated number of Ns in the assembly file | 5000 | Optional |
| mercury | **purpose_of_sequencing** | String | Populate to overwrite `purpose_of_sequencing` column |  | Optional |
| mercury | **seq_platform** | String | Populate to overwrite `seq_platform` column |  | Optional |
| mercury | **single_end** | Boolean | Set to true if your data is single-end; this ensures that a read2 column is not included in the metadata | FALSE | Optional |
| mercury | **skip_county** | Boolean | Use if your Terra table contains a county column that you do not want to include in your submission.  | FALSE | Optional |
| mercury | **state** | String | Populate to overwrite `state` column |  | Optional |
| mercury | **submitting_lab** | String | Populate to overwrite `submitting_lab` column |  | Optional |
| mercury | **submitting_lab_address** | String | Populate to overwrite `submitting_lab_address` column |  | Optional |
| mercury | **submitter_email** | String | Populate to overwrite `submitter_email` column |  | Optional |
| mercury | **usa_territory** | Boolean | If true, the "state" column will be used in place of the "country" column. For example, if "state" is Puerto Rico, then the GISAID virus name will be `hCoV-19/Puerto Rico/<name>/<year>`. The NCBI `geo_loc_name` will be "USA: Puerto Rico". This optional Boolean variable should only be used with clear understanding of what it does. | FALSE | Optional |
| mercury | **using_clearlabs_data** | Boolean | When set to `true` will change `read1_dehosted` → `clearlabs_fastq_gz`; `assembly_fasta` → `clearlabs_fasta`; `assembly_mean_coverage` → `clearlabs_sequencing_depth` | FALSE | Optional |
| mercury | **using_reads_dehosted** | Boolean | When set to true will only change read1_dehosted → reads_dehosted. Takes priority over the replacement for read1_dehosted made with the using_clearlabs_data Boolean input | FALSE | Optional |
| mercury | **vadr_alert_limit** | Int | Only for "sars-cov-2" submissions; used to filter out any samples that contain more than the indicated number of vadr alerts | 0 | Optional |
| mercury_prep_n_batch | **authors_sbt** | File | Only for "mpox" submissions; a file that contains author information. This file can be created here: <https://submit.ncbi.nlm.nih.gov/genbank/template/submission/> |  | Optional |
| mercury_prep_n_batch | **organism** | String | The organism that you want submission prepare for — each organism requires different metadata fields so please ensure this field is accurate. Options: "flu", "mpox"" or "sars-cov-2". NOTE: See **metadata_organism** for populating the `organism` column. | sars-cov-2 | Optional |
| mercury_prep_n_batch | **output_name** | String | Free text prefix for all output files | mercury | Optional |
| mercury_prep_n_batch | **skip_ncbi** | Boolean | Set to true if you only want to prepare GISAID submission files | FALSE | Optional |
| table2asn | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional |
| table2asn | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| table2asn | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/ncbi-table2asn:1.26.678 | Optional |
| table2asn | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 1 | Optional |
| trim_genbank_fastas | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional |
| trim_genbank_fastas | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| trim_genbank_fastas | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/vadr:1.3 | Optional |
| trim_genbank_fastas | **max_length** | Int | Only for "sars-cov-2" submissions; the maximum genome length for trimming terminal ambiguous nucleotides. If your sample's genome is higher than this value, the workflow will error/fail. | 30000 | Optional |
| trim_genbank_fastas | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional |
| trim_genbank_fastas | **min_length** | Int | Only for "sars-cov-2" submissions; the minimum genome length for trimming terminal ambiguous nucleotides. If your sample's genome is lower than this value, the workflow will error/fail.  | 50 | Optional |
| version_capture | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional |

</div>

### Outputs

<div class="searchable-table" markdown="1">

| **Variable** | **Type** | **Description** |
|---|---|---|
| bankit_sqn_to_email | File | **Only for mpox submission**: the sqn file that you will use to submit mpox assembly files to NCBI via email |
| biosample_metadata | File | BioSample metadata TSV file for upload to NCBI |
| excluded_samples | File | A file that contains the names and reasons why a sample was excluded from submission. **For SARS-CoV-2**, there are two sections: First, a section for any samples that failed to meet pre-determined quality thresholds (`number_N` and `vadr_num_alert`). Second, a section that includes a table that describes any missing required metadata for each sample. This table has the sample name for rows and any columns that have missing metadata as headers. If a sample is missing a piece of required metadata, the corresponding cell will be blank. However, if a different sample does have metadata for that column, the associated value will appear in the corresponding cell. **For flu and mpox**, only the second section described above exists. _Please see the example below for more details_. |
| genbank_fasta | File | **Only for SARS-CoV-2 submission**: GenBank fasta file for upload |
| genbank_metadata | File | **Only for SARS-CoV-2 submission**: GenBank metadata for upload |
| gisaid_fasta | File | **Only for mpox and SARS-CoV-2 submission**: GISAID fasta file for upload |
| gisaid_metadata | File | **Only for mpox and SARS-CoV-2 submission**: GISAID metadata for upload |
| mercury_prep_n_batch_analysis_date | String | Date analysis was run |
| mercury_prep_n_batch_version | String | Version of the PHB repository that hosts this workflow |
| mercury_script_version | String | Version of the Mercury tool that was used in this workflow |
| sra_metadata | File | SRA metadata TSV file for upload |

</div>

???+ toggle "An example excluded_samples.tsv file"

    ##### An example excluded_samples.tsv file {#example-excluded-samples}

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
