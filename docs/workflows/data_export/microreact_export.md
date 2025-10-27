---
title: Microreact_Export
---

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Name": "[**Microreact_Export**](../workflows/data_export/microreact_export.md)"}, columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level", "Dockstore"]) }}

## Microreact_Export_PHB

This set-level workflow allows for users to export Terra tables of any data directly to Microreact. If a user does not have an access token the workflow will provide a JSON that the user will be able to upload to Microreact manually. 

### Inputs

This workflow runs on the _set_ level.

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "Microreact_Export"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "Microreact_Export"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///
