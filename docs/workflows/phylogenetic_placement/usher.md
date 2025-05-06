# Usher

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Phylogenetic Placement](../../workflows_overview/workflows_type.md/#phylogenetic-placement) | [Viral](../../workflows_overview/workflows_kingdom.md/#viral) | PHB v2.1.0 | Yes | Sample-level, set-level |

## Usher_PHB

[UShER](https://usher-wiki.readthedocs.io/en/latest/) (Ultrafast Sample Placement on Existing Trees) rapidly places new samples onto an existing phylogeny using maximum parsimony. This workflow uses the UCSC-maintained global trees for SARS-CoV-2, mpox, RSV-A, and RSV-B if those organisms are specified in the `organism` input field. However, UShER can be used on any organism as long as a mutation-annotated tree (MAT) is provided in protobuf format. Contact us if you need help generating your own mutation-annotated tree, or follow the instructions available on the UShER wiki [here](https://usher-wiki.readthedocs.io/en/latest/).

### Inputs

While this workflow is technically a set-level workflow, it works on the sample-level too. When run on the set-level, the samples are placed with respect to each other.

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="Usher", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Outputs

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/tables/all_outputs.tsv", input_table=False, filter_column="Workflow", filter_values="Usher", columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///
