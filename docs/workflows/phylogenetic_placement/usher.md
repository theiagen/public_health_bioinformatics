# Usher

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Name": "[**Usher**](../workflows/phylogenetic_placement/usher.md)"}, columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level", "Dockstore"]) }}

## Usher_PHB

[UShER](https://usher-wiki.readthedocs.io/en/latest/) (Ultrafast Sample Placement on Existing Trees) rapidly places new samples onto an existing phylogeny using maximum parsimony. This workflow uses the UCSC-maintained global trees for SARS-CoV-2, mpox, RSV-A, and RSV-B if those organisms are specified in the `organism` input field. However, UShER can be used on any organism as long as a mutation-annotated tree (MAT) is provided in protobuf format. Contact us if you need help generating your own mutation-annotated tree, or follow the instructions available on the UShER wiki [here](https://usher-wiki.readthedocs.io/en/latest/).

### Inputs

While this workflow is technically a set-level workflow, it works on the sample-level too. When run on the set-level, the samples are placed with respect to each other.

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "Usher"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "Usher"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///
