---
title: Zip_Column_Content
---

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Exporting Data From Terra](../../workflows_overview/workflows_type.md/#exporting-data-from-terra) | [Any taxa](../../workflows_overview/workflows_kingdom.md/#any-taxa) | PHB v2.1.0 | Yes | Set-level |

## Zip_Column_Content_PHB

This workflow will create a zip file containing all of the items from a given column in a Terra Data Table. This is useful when you want to share a collection of result files.

### Inputs

This workflow runs on the _set_ level.

!!! tip "Do **not** include an extension"
    When you provide the name of the otutput file, **do not** include the extension. The `.zip` extension will be added to _any_ name you provide.

    While you _can_ add an extension to the `zipped_file_name` input variable, it will not affect the output file format. The output file will always be a `.zip` file.

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "Zip_Column_Content"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Outputs

!!! info "Prevent Output Overwriting"
    Please note that if you run this workflow on the _same_ Terra set, the results will overwrite each other. We recommend either (1) renaming the output variable, or (2) creating a new set every time you run the workflow. Multiple sets containing the same samples can be created as long as the set names are unique.

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "Zip_Column_Content"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///
