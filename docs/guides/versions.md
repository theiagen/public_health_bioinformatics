# Database and Docker Versions

Throughout this repository, we have carefully selected a number of default databases and Docker images that are intended to be sufficient starting points for most analyses. The following tables describe the default values for all databases, Docker images, and reference files used throughout the suite of PHB workflows.

## Default Databases

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Classification": "database"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Workflow"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

## Default Docker Images

All Docker variables are **String** types with the description: "The Docker container to use for the task".

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Classification": "docker"}, columns=["Terra Task Name", "Variable", "Default Value", "Workflow"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

## Default References

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Classification": "reference"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Workflow"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///
