# ChroQueTas 

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Name": "[**ChroQueTas**](../workflows/standalone/chroquetas.md)"}, columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level", "Dockstore"]) }}

## ChroQueTas_PHB

The ChroQueTas workflow is a standalone version of the FungAMR-integrated antimicrobial resistance screening tool. Please see the ["workflow tasks"](#workflow-tasks) section for a more detailed explanation.

### Inputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "ChroQueTas"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Workflow Tasks

{{ include_md("common_text/chroquetas_task.md") }}

### Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "ChroQueTas"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///

## References

> Bédard, C., Pageau, A., Fijarczyk, A. et al. FungAMR: a comprehensive database for investigating fungal mutations associated with antimicrobial resistance. Nat Microbiol 10, 2338–2352 (2025). https://doi.org/10.1038/s41564-025-02084-7