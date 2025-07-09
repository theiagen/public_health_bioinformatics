# TBProfiler_tNGS

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Name": "[**TBProfiler_tNGS**](../workflows/standalone/tbprofiler_tngs.md)"}, columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level"]) }}

## TBProfiler_tNGS_PHB

This workflow is still in experimental research stages. Documentation is minimal as changes may occur in the code; it will be fleshed out when a stable state has been achieved.

### Inputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "TBProfiler_tNGS"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Terra Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "TBProfiler_tNGS"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///
