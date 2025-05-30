---
title: Workflows by Kingdom
---

[Sort by Type](workflows_type.md) | [Sort Alphabetically](workflows_alphabetically.md)

---

### Any Taxa

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filter_column="Applicable Kingdom", filter_values="[Any taxa](../../workflows_overview/workflows_kingdom.md#any-taxa)", columns=["Name", "Description", "Workflow Type", "Workflow Level", "Command-line Compatibility", "Last Known Changes", "Dockstore"], replacements={"../../workflows_overview/": ""}) }}

///

### Bacteria

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filter_column="Applicable Kingdom", filter_values="[Bacteria](../../workflows_overview/workflows_kingdom.md#bacteria)", columns=["Name", "Description", "Workflow Type", "Workflow Level", "Command-line Compatibility", "Last Known Changes", "Dockstore"], replacements={"../../workflows_overview/": ""}) }}

///

### Mycotics

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filter_column="Applicable Kingdom", filter_values="[Mycotics](../../workflows_overview/workflows_kingdom.md#mycotics)", columns=["Name", "Description", "Workflow Type", "Workflow Level", "Command-line Compatibility", "Last Known Changes", "Dockstore"], replacements={"../../workflows_overview/": ""}) }}

///

### Viral

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filter_column="Applicable Kingdom", filter_values="[Viral](../../workflows_overview/workflows_kingdom.md#viral)", columns=["Name", "Description", "Workflow Type", "Workflow Level", "Command-line Compatibility", "Last Known Changes", "Dockstore"], replacements={"../../workflows_overview/": ""}) }}

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
