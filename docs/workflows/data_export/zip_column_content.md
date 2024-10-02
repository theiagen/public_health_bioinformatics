---
title: Zip_Column_Content
---

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Exporting Data From Terra](../../workflows_overview/workflows_type.md/#exporting-data-from-terra) | [Any taxa](../../workflows_overview/workflows_kingdom.md/#any-taxa) | PHB v2.1.0 | Yes | Set-level |

## Zip_Column_Content_PHB

This workflow will create a zip file that contains all of the items in a column in a Terra Table.

### Inputs

This workflow runs on the set level.

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| zip_column_content | **files_to_zip** | Array[File] | The column that has the files you want to zip. |  | Required |
| zip_column_content | **zipped_file_name** | String | The name you want your zipped file to have. The .zip file extension will be added to this name. |  | Required |
| zip_files | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| zip_files | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| zip_files | **docker_image** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1" | Optional |
| zip_files | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| version_capture | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional |

### Outputs

!!! info
    Please note that if you run this workflow on the same Terra set (the same group of samples can be included in multiple Terra sets), the results will overwrite each other. We recommend either (1) renaming the output variable, or (2) creating a new set every time you run the workflow.

| **Variable** | **Type** | **Description** |
|---|---|---|
| zipped_files | File | The zipped file containing all of the items from the column you selected. |
| zip_column_content_version | String | The version of the repository the workflow is hosted in |
| zip_column_content_analysis_date | String | The date the workflow was run |
