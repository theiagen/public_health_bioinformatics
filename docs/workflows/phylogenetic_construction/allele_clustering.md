# Workflow Name

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Name": "[**Allele Clustering**](../workflows/phylogenetic_construction/allele_clustering.md)"}, columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level", "Dockstore"]) }}

## Allele_Clustering_PHB

The Allele Clustering module is used by PulseNet 2.0 to generate NWK trees for visualization, using the results from the `allele_clustering` task in TheiaProk.

### Inputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "Allele_Clustering"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Workflow Tasks

{{ include_md("common_text/versioning_task.md") }}
{{ include_md("common_text/allele_clustering_task.md") }}

### Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "Allele_Clustering"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

//

## References

> [pulsenet2.0-trees](https://github.com/ncezid-biome/pulsenet2.0-trees/tree/main)
