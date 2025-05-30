---
title: Concatenate_Column_Content
---

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Exporting Data From Terra](../../workflows_overview/workflows_type.md/#exporting-data-from-terra) | [Any taxa](../../workflows_overview/workflows_kingdom.md/#any-taxa) | PHB v2.1.0 | Yes | Set-level |

## Concatenate_Column_Content_PHB

This set-level workflow will create a file containing all of the items from a given column in a Terra Data Table. This is useful when you want to investigate a collection of result files. There is a video available with more information about the Concatenate_Column_Content workflow: **📺 [Workflow Focus: Concatenate_Column_Content](https://www.youtube.com/watch?v=T5Gnj9BtC9I)**. Although this video refers to an older version of this workflow and various names may be changed, the concepts presented are still applicable.

### Inputs

This workflow runs on the _set_ level.

!!! tip "Include the extension"
    When you provide the name of the output file, be sure to include the extension. For example, if you want the output file to be a FASTA file, you should include the `.fasta` extension in the `concatenated_file_name` input variable. Otherwise, the workflow will create a file without any extension, which can cause problems when you try to open it with certain programs and operating systems.

    Please note that this workflow **cannot** produce Excel files. If you need an Excel file, you can convert the output file to Excel using a program like Microsoft Excel or Google Sheets.

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "Concatenate_Column_Content"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Outputs

!!! info "Prevent Output Overwriting"
    Please note that if you run this workflow on the _same_ Terra set, the results will overwrite each other. We recommend either (1) renaming the output variable, or (2) creating a new set every time you run the workflow. Multiple sets containing the same samples can be created as long as the set names are unique.

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "Concatenate_Column_Content"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///
