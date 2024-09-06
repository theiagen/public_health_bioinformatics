---
title: Transfer_Column_Content
---

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** ||
|---|---|---|---|---|
| [Exporting Data From Terra](../../workflows_overview/workflows-type.md/#exporting-data-from-terra) | [Any taxa](../../workflows_overview/workflows-kingdom.md/#any-taxa) | PHB v1.3.0 | Yes | Set-level |

## Transfer_Column_Content_PHB

This set-level workflow will transfer all of the items from a given column in a Terra Data Table to a single GCP storage bucket location. This is useful when you want to transfer many files to another GCP storage bucket (can be a Terra workspace storage bucket or a non-Terra storage bucket). 

!!! note
    This workflow requires that the user's Terra pet-service account has sufficient privileges to read and write to the target storage bucket.

    - If the target bucket **is associated with a Terra workspace**, the workspace OWNER/administrator must grant WRITER privileges with the Terra workspace.
    - If the target bucket **is not associated with a Terra workspace** (i.e. GCP storage bucket), the user's Terra pet-service account (or their Terra PROXY account) must be granted the ability to read and write to the bucket (Storage Object Admin google privileges)

!!! note
    If using Transfer_column_content workflow version 1.3.0 or higher, the call-caching feature of Terra has been DISABLED to ensure that the workflow is run from the beginning and data is transferred fresh. Call-caching will not be enabled, even if the user checks the box ✅ in the Terra workflow interface.

### Inputs

This workflow runs on the set level.

| **Terra Task name** | **input_variable** | **Type** | **Description** | **Default attribute** | **Status** |
|---|---|---|---|---|---|
| transfer_column_content | **files_to_transfer** | Array[File] | The column that has the files you want to concatenate. | | Required |
| transfer_column_content | **target_bucket** | String | The GS URI of the target storage bucket. Note: **Do not include spaces**, but **do** include the `gs://` at the beginning of the bucket URI | | Required |
| transfer_files | **cpu** | Int | Number of cpus to allocate to this task | 4 | Optional |
| transfer_files | **disk_size** | Int | Storage in GB for the disk for this task | 100 | Optional |
| transfer_files | **docker_image** | String | The docker image used to perform the file transfer. | us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1" | Optional |
| transfer_files | **memory** | Int | RAM in GB for this task | 8 | Optional |
| version_capture | **docker** | String | The Docker image used to run the version_capture task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional |

### Outputs

!!! info
    Please note that if you run this workflow on the same Terra set (the same group of samples can be included in multiple Terra sets), the results will overwrite each other. We recommend either (1) renaming the output variable, or (2) creating a new set every time you run the workflow.

| **Variable** | **Type** | **Description** |
|---|---|---|
| transferred_files | File | A list of all of the files now located at the target bucket location (GSURI) |
| transfer_column_content_version | String | The version of the repository the workflow is hosted in |
| transfer_column_content_analysis_date | String | The date the workflow was run |
