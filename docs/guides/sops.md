# Available Standard Operating Procedures

Theiagen has developed a number of Standard Operating Procedures (SOPs) to help you analyze your data. Please see the available SOPs below.

## All SOPs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_sops.tsv", sort_by=["Current or prior wf version?", "SOP"], columns=["SOP", "SOP version", "Workflow", "PHB version compatibility", "Pathogen/Category", "Current or prior wf version?"]) }}

///

## Most Up-to-Date SOPs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_sops.tsv", sort_by=["SOP", "SOP version"], filters={"Current or prior wf version?": "Current"}, columns=["SOP", "SOP version", "Workflow", "PHB version compatibility", "Pathogen/Category", "Current or prior wf version?"]) }}

///

## Data Import SOPs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_sops.tsv", sort_by=["Current or prior wf version?", "SOP"], filters={"Pathogen/Category": "Data Import"}, columns=["SOP", "SOP version", "Workflow", "PHB version compatibility", "Pathogen/Category", "Current or prior wf version?"]) }}

///

## Viral SOPs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_sops.tsv", sort_by=["Current or prior wf version?", "SOP"], filters={"Pathogen/Category": ["Flu", "SC2"]}, columns=["SOP", "SOP version", "Workflow", "PHB version compatibility", "Pathogen/Category", "Current or prior wf version?"]) }}

///

## Bacterial SOPs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_sops.tsv", sort_by=["Current or prior wf version?", "SOP"], filters={"Pathogen/Category": ["Bacterial", "HAI"]}, columns=["SOP", "SOP version", "Workflow", "PHB version compatibility", "Pathogen/Category", "Current or prior wf version?"]) }}

///

## Phylogenetic SOPs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_sops.tsv", sort_by=["Current or prior wf version?", "SOP"], filters={"Pathogen/Category": "Phylogenetic construction"}, columns=["SOP", "SOP version", "Workflow", "PHB version compatibility", "Pathogen/Category", "Current or prior wf version?"]) }}

///

## Data Sharing SOPs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_sops.tsv", sort_by=["Current or prior wf version?", "SOP"], filters={"Pathogen/Category": "Public data sharing"}, columns=["SOP", "SOP version", "Workflow", "PHB version compatibility", "Pathogen/Category", "Current or prior wf version?"]) }}

///

## Miscellaneous SOPs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_sops.tsv", sort_by=["Current or prior wf version?", "SOP"], filters={"Pathogen/Category": ["Data Import", "Downsampling", "Getting Started", "Validations"]}, columns=["SOP", "SOP version", "Workflow", "PHB version compatibility", "Pathogen/Category", "Current or prior wf version?"]) }}

///
