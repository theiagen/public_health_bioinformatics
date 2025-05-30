---
title: Transfer_Column_Content
---

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filter_column="Name", filter_values="[**Transfer_Column_Content**](../workflows/data_export/transfer_column_content.md)", columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level"]) }}

## Transfer_Column_Content_PHB

This set-level workflow will transfer all of the files from a given column in a Terra Data Table to a single Google Cloud Platform (GCP) storage bucket location. This is useful when you want to transfer many files to another GCP storage bucket (this can be a Terra workspace storage bucket or a non-Terra storage bucket).

!!! tip "Ensure Proper Permissions"
    This workflow requires that the user's Terra pet-service account has sufficient privileges to read and write to the target storage bucket.

    - If the target bucket **is associated with a Terra workspace**, the workspace OWNER/administrator must grant you WRITER privileges to the Terra workspace.
    - If the target bucket **is not associated with a Terra workspace** (i.e. GCP storage bucket), the user's Terra pet-service account (or their Terra PROXY account) must be granted the ability to read and write to the bucket (Storage Object Admin Google privileges)

!!! info "Call-Caching Disabled"
    If using Transfer_Column_Content workflow version 1.3.0 or higher, the call-caching feature of Terra has been DISABLED to ensure that the workflow is run from the beginning and data is transferred fresh. Call-caching will not be enabled, even if the user checks the box ✅ in the Terra workflow interface.

### Inputs

This workflow runs on the _set_ level.

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "Transfer_Column_Content"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

[^1]: GS URI: Google Storage Uniform Resource Identifier. This is **not** the same as a URL, which typically begins with http:// or https://. A GS URI begins with `gs://` and is used to reference a location in a Google Cloud Storage Bucket. For example, `gs://bucket-name/folder-name/file-name`. Other cloud storage providers have their own URIs, such as `s3://` for Amazon S3, although this workflow only supports Google Cloud Storage URIs.

### Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "Transfer_Column_Content"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///
