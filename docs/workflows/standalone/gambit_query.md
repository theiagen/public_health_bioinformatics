# GAMBIT_Query

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filter_column="Name", filter_values="[**GAMBIT_Query**](../workflows/standalone/gambit_query.md)", columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level"]) }}

## GAMBIT_Query_PHB

The GAMBIT_Query_PHB workflow performs taxon assignment of a genome assembly using the GAMBIT task.

### Inputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="Gambit_Query", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Workflow Tasks

{{ include_md("common_text/gambit_task.md") }}

### Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filter_column="Workflow", filter_values="Gambit_Query", columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///

> GAMBIT (Genomic Approximation Method for Bacterial Identification and Tracking): A methodology to rapidly leverage whole genome sequencing of bacterial isolates for clinical identification. Lumpe et al. PLOS ONE, 2022. DOI: [10.1371/journal.pone.0277575](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0277575)
