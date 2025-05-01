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

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/input_tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="CZGenEpi_Prep", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

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
