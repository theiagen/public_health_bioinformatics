---
title: Workflows by Kingdom
---

[Sort by Type](workflows_type.md) | [Sort Alphabetically](workflows_alphabetically.md)

---

### Any Taxa

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filter_column="Applicable Kingdom", filter_values="Any taxa", columns=["Name", "Description", "Applicable Kingdom", "Workflow Level", "Command-line Compatibility", "Last Known Changes", "Dockstore"]) }}

///

### Bacteria

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filter_column="Applicable Kingdom", filter_values="Bacteria", columns=["Name", "Description", "Applicable Kingdom", "Workflow Level", "Command-line Compatibility", "Last Known Changes", "Dockstore"]) }}

///

### Mycotics

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filter_column="Applicable Kingdom", filter_values="Mycotics", columns=["Name", "Description", "Applicable Kingdom", "Workflow Level", "Command-line Compatibility", "Last Known Changes", "Dockstore"]) }}

///

### Viral

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filter_column="Applicable Kingdom", filter_values="Viral", columns=["Name", "Description", "Applicable Kingdom", "Workflow Level", "Command-line Compatibility", "Last Known Changes", "Dockstore"]) }}

///

<!-- definitions for workflow type column -->
*[Sample-level]: This workflow is run once for each sample
*[Set-level]: This workflow is run once on a group of samples

<!-- definitions for taxa column -->
*[Any taxa]: This workflow is organism-agnostic and can be run with any taxa
*[Viral]: This workflow is compatible with any viral pathogen
*[Bacteria]: This workflow is compatible with any bacterial pathogen
*[Mycotics]: This workflow is compatible with mycotic pathogens

<!-- definition for command-line compatible column -->
*[Command-line Compatibility]: Command-line compatibility is determined if the workflow can be run on a local command-line environment, providing all dependencies are installed, with either `miniwdl` or `cromwell`.
*[Some optional features incompatible]: Some optional features of this workflow are incompatible with command-line use and require modification
*[Yes]: This workflow is compatible with command-line use
*[No]: This workflow is not compatible with command-line use even with modifications
