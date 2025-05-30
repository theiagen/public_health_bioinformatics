# Terra_2_GISAID

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filter_column="Name", filter_values="[**Terra_2_GISAID**](../workflows/public_data_sharing/terra_2_gisaid.md)", columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level"]) }}

## Terra_2_GISAID_PHB

Terra_2_GISAID programmatically submits SARS-CoV-2 assembly files to GISAID.

This workflow expects data that has been prepared for submission using either Mercury_Batch or Mercury_Prep_N_Batch (recommended).

!!! dna "client-ID"
    To obtain a client-ID, contact `clisupport@gisaid.org` and include your username in your request.

### Inputs

The optional variable `frameshift_notification` has three options that correspond to the associated web-browser options:

- "**catch_all**" - "Notify me about ALL DETECTED FRAMESHIFTS in this submission for reconfirmation of affected sequences"
- "**catch_novel**" [DEFAULT] - "Notify me only about NOT PREVIOUSLY REPORTED FRAMESHIFTS in this submission for reconfirmation of affected sequences"
- "**catch_none**" - "I confirm ANY FRAMESHIFTS in this submission and request their release without confirmation by a curator"

!!! warning "GISAID Credentials"
    Please note that the user must provide either an authentication_file or a gisaid_credentials file to run this workflow; explanations for both can be found in the table below.

This workflow runs on the sample level.

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="Terra_2_GISAID", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filter_column="Workflow", filter_values="Terra_2_GISAID", columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///
