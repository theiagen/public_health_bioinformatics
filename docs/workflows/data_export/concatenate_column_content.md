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

<div class="searchable-table" markdown="1">

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| concatenate_column_content | **concatenated_file_name** | String | The name of the output file. **_Include the extension_**, such as ".fasta" or ".txt". |  | Required |
| concatenate_column_content | **files_to_cat** | Array[File] | The column that has the files you want to concatenate. |  | Required |
| cat_files | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| cat_files | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| cat_files | **docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1 | Optional |
| cat_files | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| cat_files | **skip_extra_headers** | Boolean | If the files you are concatenating have identical headers, you can include only the first instance of the header and skip all of the others so they do not appear duplicated in the concatenated file. To activate this, set to true. | false | Optional |
| version_capture | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0 | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional |

</div>

### Outputs

!!! info "Prevent Output Overwriting"
    Please note that if you run this workflow on the _same_ Terra set, the results will overwrite each other. We recommend either (1) renaming the output variable, or (2) creating a new set every time you run the workflow. Multiple sets containing the same samples can be created as long as the set names are unique.

| **Variable** | **Type** | **Description** |
|---|---|---|
| concatenated_files | File | The file containing all of the items from the column you selected. |
| concatenate_column_content_version | String | The version of the repository the workflow is hosted in |
| concatenate_column_content_analysis_date | String | The date the workflow was run |
