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

<div class="searchable-table" markdown="1">

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| zip_column_content | **files_to_zip** | Array[File] | The column that has the files you want to zip. |  | Required |
| zip_column_content | **zipped_file_name** | String | The name you want your zipped file to have. The .zip file extension will be added to this name. |  | Required |
| zip_files | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| zip_files | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| zip_files | **docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1 | Optional |
| zip_files | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| version_capture | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0 | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional |

</div>

### Outputs

!!! info "Prevent Output Overwriting"
    Please note that if you run this workflow on the _same_ Terra set, the results will overwrite each other. We recommend either (1) renaming the output variable, or (2) creating a new set every time you run the workflow. Multiple sets containing the same samples can be created as long as the set names are unique.

| **Variable** | **Type** | **Description** |
|---|---|---|
| zipped_files | File | The zipped file containing all of the items from the column you selected. |
| zip_column_content_version | String | The version of the repository the workflow is hosted in |
| zip_column_content_analysis_date | String | The date the workflow was run |
