---
title: Workflows by Kingdom
---

[Sort by Type](workflows_type.md) | [Sort Alphabetically](workflows_alphabetically.md)

---

### Any taxa

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Applicable Kingdom": "[Any taxa](../../workflows_overview/workflows_kingdom.md#any-taxa)"}, columns=["Name", "Description", "Workflow Type", "Workflow Level", "Command-line Compatibility", "Last Known Changes", "Dockstore"], replacements={"../../workflows_overview/": ""}) }}

///

### Bacteria

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Applicable Kingdom": "[Bacteria](../../workflows_overview/workflows_kingdom.md#bacteria)"}, columns=["Name", "Description", "Workflow Type", "Workflow Level", "Command-line Compatibility", "Last Known Changes", "Dockstore"], replacements={"../../workflows_overview/": ""}) }}

///

### Mycotics

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Applicable Kingdom": "[Mycotics](../../workflows_overview/workflows_kingdom.md#mycotics)"}, columns=["Name", "Description", "Workflow Type", "Workflow Level", "Command-line Compatibility", "Last Known Changes", "Dockstore"], replacements={"../../workflows_overview/": ""}) }}

///

### Viral

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Applicable Kingdom": "[Viral](../../workflows_overview/workflows_kingdom.md#viral)"}, columns=["Name", "Description", "Workflow Type", "Workflow Level", "Command-line Compatibility", "Last Known Changes", "Dockstore"], replacements={"../../workflows_overview/": ""}) }}

///

<!-- definitions for workflow level column -->
*[Sample-level]: These workflows are run once for each sample
*[Set-level]: These workflows are run once on a group of samples
*[Table-level]: These workflows create or process entire tables of data.

<!-- definitions for applicable kingdom column -->
*[Any taxa]: These workflows are organism-agnostic and can be run with any taxa
*[Viral]: These workflows are compatible with any viral pathogen
*[Bacteria]: These workflows are compatible with any bacterial pathogen
*[Mycotics]: These workflows are compatible with mycotic pathogens

<!-- definitions for workflow type column -->
*[Comparative Analysis]: These workflows compare results between analyses
*[Data Import]: These workflows either import data to Terra or perform file manipulation prior to analysis.
*[Genomic Characterization]: These workflows perform genomic characterization of samples.
*[Phylogenetic Construction]: These workflows construct phylogenetic trees or assist in that process.
*[Phylogenetic Placement]: These workflows place samples on existing phylogenetic trees.
*[Public Data Sharing]: These workflows help prepare and/or transfer data to public data repositories.
*[Exporting Data from Terra]: These workflows export data from Terra for use in other platforms or tools.
*[Standalone]: These workflows are designed to be run independently of the major workflow groups as either supplements or alternatives.

<!-- definition for command-line compatible column -->
*[Command-line Compatibility]: Command-line compatibility is determined if the workflow can be run on a local command-line environment, providing all dependencies are installed, with either `miniwdl` or `cromwell`.
*[Some optional features incompatible]: Some optional features of these workflows are incompatible with command-line use and require modification
*[Yes]: These workflows are compatible with command-line use
*[No]: These workflows are not compatible with command-line use even with modifications; these workflows are engineered to run on Terra
