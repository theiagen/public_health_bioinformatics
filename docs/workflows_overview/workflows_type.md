---
title: Workflows by Type
---

[Sort by Kingdom](workflows_kingdom.md) | [Sort Alphabetically](workflows_alphabetically.md)

---

### Data Import

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filter_column="Workflow Type", filter_values="[Data Import](../../workflows_overview/workflows_type.md#data-import)", columns=["Name", "Description", "Applicable Kingdom", "Workflow Level", "Command-line Compatibility", "Last Known Changes", "Dockstore"], replacements={"../../workflows_overview/": ""}) }}

///

### Genomic Characterization

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filter_column="Workflow Type", filter_values="[Genomic Characterization](../../workflows_overview/workflows_type.md#genomic-characterization)", columns=["Name", "Description", "Applicable Kingdom", "Workflow Level", "Command-line Compatibility", "Last Known Changes", "Dockstore"], replacements={"../../workflows_overview/": ""}) }}

///

### Phylogenetic Construction

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filter_column="Workflow Type", filter_values="[Phylogenetic Construction](../../workflows_overview/workflows_type.md#phylogenetic-construction)", columns=["Name", "Description", "Applicable Kingdom", "Workflow Level", "Command-line Compatibility", "Last Known Changes", "Dockstore"], replacements={"../../workflows_overview/": ""}) }}

///

### Phylogenetic Placement

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filter_column="Workflow Type", filter_values="[Phylogenetic Placement](../../workflows_overview/workflows_type.md#phylogenetic-placement)", columns=["Name", "Description", "Applicable Kingdom", "Workflow Level", "Command-line Compatibility", "Last Known Changes", "Dockstore"], replacements={"../../workflows_overview/": ""}) }}

///

### Public Data Sharing

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filter_column="Workflow Type", filter_values="[Public Data Sharing](../../workflows_overview/workflows_type.md#public-data-sharing)", columns=["Name", "Description", "Applicable Kingdom", "Workflow Level", "Command-line Compatibility", "Last Known Changes", "Dockstore"], replacements={"../../workflows_overview/": ""}) }}

///

### Exporting Data from Terra

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filter_column="Workflow Type", filter_values="[Exporting Data from Terra](../../workflows_overview/workflows_type.md#exporting-data-from-terra)", columns=["Name", "Description", "Applicable Kingdom", "Workflow Level", "Command-line Compatibility", "Last Known Changes", "Dockstore"], replacements={"../../workflows_overview/": ""}) }}

///

### Standalone

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filter_column="Workflow Type", filter_values="[Standalone](../../workflows_overview/workflows_type.md#standalone)", columns=["Name", "Description", "Applicable Kingdom", "Workflow Level", "Command-line Compatibility", "Last Known Changes", "Dockstore"], replacements={"../../workflows_overview/": ""}) }}

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
