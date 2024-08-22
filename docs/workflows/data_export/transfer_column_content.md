---
title: Transfer_Column_Content
---

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line compatibliity** | **Workflow type** |
|---|---|---|---|---|
| [Exporting Data From Terra](../../../workflows_overview/workflows-type/#exporting-data-from-terra) | [Any taxa](../../../workflows_overview/workflows-kingdom/#any-taxa) | PHB v1.3.0 | Yes | Set-level |

## Transfer_Column_Content_PHB

This set-level workflow will transfer all of the items from a given column in a Terra Data Table to a single GCP storage bucket location. This is useful when you want to transfer many files to another GCP storage bucket (can be a Terra workspace storage bucket or a non-Terra storage bucket). 

!!! note
    This workflow requires that the user’s Terra pet-service account has sufficient privileges to read and write to the target storage bucket.

    - If the target bucket **is associated with a Terra workspace**, the workspace OWNER/administrator must grant WRITER privileges with the Terra workspace.
    - If the target bucket **is not associated with a Terra workspace** (i.e. GCP storage bucket), the user’s Terra pet-service account (or their Terra PROXY account) must be granted the ability to read and write to the bucket (Storage Object Admin google privileges)

!!! note
    If using Transfer_column_content workflow version 1.3.0 or higher, the call-caching feature of Terra has been DISABLED to ensure that the workflow is run from the beginning and data is transferred fresh. Call-caching will not be enabled, even if the user checks the box ✅ in the Terra workflow interface.

## Inputs

| **Terra Task name** | **input_variable** | **Type** | **Description** | **Status** |
|---|---|---|---|---|
| transfer_column_content | **files_to_transfer** | Array[File] | The column that has the files you want to concatenate. | Required |
| transfer_column_content | **target_bucket** | String | The GS URI of the target storage bucket. <br>Note: DO NOT INCLUDE SPACES. <br>Note: DO INCLUDE THE gs:// at the beginning | Required |
| transfer_files | **cpu** | Int | Number of cpus to allocate to this task | Optional  |
| transfer_files | **disk_size** | Int | Storage in GB for the disk for this task | Optional  |
| transfer_files | **docker_image** | String | The docker image used to perform the file transfer. | Optional  |
| transfer_files | **memory** | Int | RAM in GB for this task | Optional  |
| version_capture | **docker** | String | Docker image used to run version_capture task | Optional  |
| version_capture | **timezone** | String | The timezone  | Optional |
