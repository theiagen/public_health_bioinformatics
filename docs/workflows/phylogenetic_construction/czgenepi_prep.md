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

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="CZGenEpi_Prep", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Outputs

The concatenated_czgenepi_fasta and concatenated_czgenepi_metadata files can be uploaded directly to CZ GEN EPI without any adjustments.

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filter_column="Workflow", filter_values="CZGenEpi_Prep", columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///

## References

> CZ GEN EPI Help Center "Uploading Data" <https://help.czgenepi.org/hc/en-us/articles/6160372401172-Uploading-data>
