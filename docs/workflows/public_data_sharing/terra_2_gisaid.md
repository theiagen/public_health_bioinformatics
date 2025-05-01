# Terra_2_GISAID

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Public Data Sharing](../../workflows_overview/workflows_type.md/#public-data-sharing) | [Viral](../../workflows_overview/workflows_kingdom.md/#viral) | PHB v1.2.1 | Yes | Set-level |

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

{{ input_table("docs/assets/input_tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="Terra_2_GISAID", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Outputs

<div class="searchable-table" markdown="1">

| **Variable** | **Type** | **Description** |
|---|---|---|
| failed_uploads | Boolean | The metadata for any failed uploads |
| gisaid_cli_version | String | The verison of the GISAID CLI tool |
| gisaid_logs | File | The log files regarding the submission |
| terra_2_gisaid_analysis_date | String | The date of the analysis |
| terra_2_gisaid_version | String | The version of the PHB repository that this workflow is hosted in |

</div>
