---
title: Concatenate_Column_Content
---

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Exporting Data From Terra](../../workflows_overview/workflows-type.md/#exporting-data-from-terra) | [Any taxa](../../workflows_overview/workflows-kingdom.md/#any-taxa) | PHB v2.1.0 | Yes | Set-level |

## Concatenate_Column_Content_PHB

This set-level workflow will create a file containing all of the items from a given column in a Terra Data Table. This is useful when you want to investigate many results files. There is a video available with more information about the Concatenate_Column_Content workflow: **📺 [Workflow Focus: Concatenate_Column_Content](https://www.youtube.com/watch?v=T5Gnj9BtC9I)**

### Inputs

This workflow runs on the set level.

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| concatenate_column_content | **concatenated_file_name** | String | The name of the output file. ***Include the extension***, such as ".fasta" or ".txt". |  | Required |
| concatenate_column_content | **files_to_cat** | Array[File] | The column that has the files you want to concatenate. |  | Required |
| cat_files | **cpu** | Int | number of cpus used to run cat_files | 2 | Optional |
| cat_files | **disk_size** | Int | Amount of storage in GigaBytes (GB) requested for the VM to run the cat_files task | 100 | Optional |
| cat_files | **docker_image** | String | Docker image used to run the cat_files task | "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1" | Optional |
| cat_files | **memory** | Int | Amount of RAM/memory requested for running the cat_files task | 8 | Optional |
| cat_files | **skip_extra_headers** | Boolean | If the files you are concatenating have identical headers, you can include only the first instance of the header and skip all of the others so they do not appear duplicated in the concatenated file. To activate this, set to true. | false | Optional |
| version_capture | **docker** | String | The Docker image used to run the version_capture task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional |

### Outputs

!!! info
    Please note that if you run this workflow on the same Terra set (the same group of samples can be included in multiple Terra sets), the results will overwrite each other. We recommend either (1) renaming the output variable, or (2) creating a new set every time you run the workflow.

| **Variable** | **Type** | **Description** |
|---|---|---|
| concatenated_files | File | The file containing all of the items from the column you selected. |
| concatenate_column_content_version | String | The version of the repository the workflow is hosted in |
| concatenate_column_content_analysis_date | String | The date the workflow was run |
