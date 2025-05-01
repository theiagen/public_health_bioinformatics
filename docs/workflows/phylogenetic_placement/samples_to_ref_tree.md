# Samples_to_Ref_Tree

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Phylogenetic Placement](../../workflows_overview/workflows_type.md/#phylogenetic-placement) | [Viral](../../workflows_overview/workflows_kingdom.md/#viral) | PHB v2.1.0 | Yes | Sample-level, Set-level |

## Samples_to_Ref_Tree_PHB

[Nextclade](https://docs.nextstrain.org/projects/nextclade/en/stable/index.html) rapidly places new samples onto an existing reference phylogenetic tree. Phylogenetic placement is done by comparing the mutations of the query sequence (relative to the reference) with the mutations of every node and tip in the reference tree, and finding the node which has the most similar set of mutations. This operation is repeated for each query sequence, until all of them are placed onto the tree. This workflow uses the Nextstrain-maintained [nextclade datasets](https://github.com/nextstrain/nextclade_data) for SARS-CoV-2, mpox, influenza A and B, and RSV-A and RSV-B. The organism must be specified as input in the field `organism`, and these align with the nextclade dataset names, i.e. " sars-cov-2", "flu_h1n1pdm_ha", "flu_h1n1pdm_na", "flu_h3n2_ha", "flu_h3n2_na", "flu_vic_ha", "flu_vic_na", "flu_yam_ha", "hMPXV", "hMPXV_B1", "MPXV", "rsv_a" and "rsv_b".

However, nextclade can be used on any organism as long as an an existing, high-quality input reference tree with representative samples on it is provided, in addition to other optional inputs. Contact us if you need help generating your own mutation-annotated tree, or follow the instructions available on the Augur wiki [here](https://docs.nextstrain.org/projects/augur/en/stable/index.html).

!!! info "_Placement_ not _construction_"
    This workflow is not for building a tree from scratch, but rather for the placement of new sequences onto an existing high-quality input reference tree with representative samples on it. In effect, query samples are only compared to reference samples and never to the other query samples.

### Inputs

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/input_tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="Samples_to_Ref_Tree", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///
### Outputs

<div class="searchable-table" markdown="1">

| **Variable** | **Type** | **Description** |
|---|---|---|
| treeUpdate_auspice_json | File | Phylogenetic tree with user placed samples |
| treeUpdate_nextclade_docker | String | Nextclade docker image used |
| treeUpdate_nextclade_json | File | JSON file with the results of the Nextclade analysis |
| treeUpdate_nextclade_tsv | File | Tab-delimited file with Nextclade results |
| treeUpdate_nextclade_version | String | Nextclade version used |
| samples_to_ref_tree_analysis_date | String | Date of analysis |
| samples_to_ref_tree_version | String | Version of the Public Health Bioinformatics (PHB) repository used |

</div>
