# CZGenEpi_Prep

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Phylogenetic Construction](../../workflows_overview/workflows_type.md/#phylogenetic-construction) | [Viral](../../workflows_overview/workflows_kingdom.md/#viral) | PHB v1.3.0 | No | Set-level |

## CZGenEpi_Prep_PHB

The CZGenEpi_Prep workflow prepares data for upload to the Chan Zuckerberg GEN EPI platform, where phylogenetic trees and additional data processing can occur. This workflow extracts the necessary metadata fields from your Terra table.

### Inputs

In order to enable customization for where certain fields should be pulled from the Terra table, the user can specify different column names in the appropriate location. For example, if the user wants to use the "clearlabs_fasta" column for the assembly file _instead_ of the default "assembly_fasta" column, they can write "clearlabs_fasta" for the `assembly_fasta_column_name` optional variable.

Variables with both the "Optional" and "Required" tag require the column (regardless of name) to be present in the data table.

This workflow runs on the set level.

<div class="searchable-table" markdown="1">

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| czgenepi_prep | **sample_names** | Array[String] | The array of sample ids you want to prepare for CZ GEN EPI |  | Required |
| czgenepi_prep | **terra_table_name** | String | The name of the Terra table where the data is hosted |  | Required |
| czgenepi_prep | **terra_project_name** | String | The name of the Terra project where the data is hosted |  | Required |
| czgenepi_prep | **terra_workspace_name** | String | The name of the Terra workspace where the data is hosted |  | Required |
| download_terra_table | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 10 | Optional |
| download_terra_table | **docker** | String | The Docker container to use for the task | quay.io/theiagen/terra-tools:2023-06-21 | Optional |
| download_terra_table | **disk_size** | String | The size of the disk used when running this task | 1 | Optional |
| download_terra_table | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional |
| czgenepi_prep | **assembly_fasta_column_name** | String | The column name where the sample's assembly file can be found | assembly_fasta | Optional, Required |
| czgenepi_prep | **county_column_name** | String | The column name where the samples' originating county can be found | county | Optional, Required |
| czgenepi_prep | **organism** | String | The organism for data preparation. Options: "mpox" or "sars-cov-2" | sars-cov-2 | Optional |
| czgenepi_prep | **is_private** | Boolean | Sets whether sample status is provate, or not | true | Optional |
| czgenepi_prep | **genbank_accession_column_name** | String | The column name where the genbank accession for the sample can be found | genbank_accession | Optional |
| czgenepi_prep | **country_column_name** | String | The column name where the sample's originating country can be found | country | Optional, Required |
| czgenepi_prep | **collection_date_column_name** | String | The column name where the sample's collection date can be found | collection_date | Optional, Required |
| czgenepi_prep | **state_column_name** | String | The column name where the sample's originating state can be found | state | Optional, Required |
| czgenepi_prep | **continent_column_name** | String | The column name where the sample's originating continent can be found | continent | Optional, Required |
| czgenepi_prep | **sequencing_date_column_name** | String | The column name where the sample's sequencing data can be found | sequencing_date | Optional |
| czgenepi_prep | **private_id_column_name** | String | The column name where the Private ID for the sample can be found | terra_table_name_id | Optional, Required |
| czgenepi_wrangling | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| czgenepi_wrangling | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-08-08-2 | Optional |
| czgenepi_wrangling | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| czgenepi_wrangling | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional |
| version_capture | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional |

</div>

### Outputs

The concatenated_czgenepi_fasta and concatenated_czgenepi_metadata files can be uploaded directly to CZ GEN EPI without any adjustments.

| **Variable** | **Type** | **Description** |
|---|---|---|
| concatenate_czgenepi_fasta | File | The concatenated fasta file with the renamed headers (the headers are renamed to account for clearlabs data which has unique headers) |
| concatenate_czgenepi_metadata | File | The concatenated metadata that was extracted from the terra table using the specified columns |
| czgenepi_prep_version | String | The version of PHB the workflow is in |
| czgenepi_prep_analysis_date | String | The date the workflow was run |

## References

> CZ GEN EPI Help Center "Uploading Data" <https://help.czgenepi.org/hc/en-us/articles/6160372401172-Uploading-data>
