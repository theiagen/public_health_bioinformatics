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

{{ input_table("docs/assets/input_tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="Usher", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Outputs

<div class="searchable-table" markdown="1">

| **Variable** | **Type** | **Description** |
|---|---|---|
| usher_clades | File | The clades predicted for the samples |
| usher_phb_analysis_date | String | The date the analysis was run |
| usher_phb_version | String | The version of PHB the workflow is from |
| usher_protobuf_version | String | The version of the mutation-annotated protobuf tree (what day and what samples are included, if a default organism was used; otherwise, says it was user-provided) |
| usher_subtree_mutations | Array[File] | An array of files showing the mutations at each internal node for the subtree |
| usher_subtrees | Array[File] | An array of subtrees where your samples have been placed |
| usher_uncondensed_tree | File | The entire global tree with your samples included (warning: may be a very large file if the organism is "sars-cov-2") |
| usher_version | String | The version of UShER used |

</div>
