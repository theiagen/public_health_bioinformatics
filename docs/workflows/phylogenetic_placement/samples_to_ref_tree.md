# Samples_to_Ref_Tree

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Name": "[**Samples_to_Ref_Tree**](../workflows/phylogenetic_placement/samples_to_ref_tree.md)"}, columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level"]) }}

## Samples_to_Ref_Tree_PHB

[Nextclade](https://docs.nextstrain.org/projects/nextclade/en/stable/index.html) rapidly places new samples onto an existing reference phylogenetic tree. Phylogenetic placement is done by comparing the mutations of the query sequence (relative to the reference) with the mutations of every node and tip in the reference tree, and finding the node which has the most similar set of mutations. This operation is repeated for each query sequence, until all of them are placed onto the tree. This workflow uses the Nextstrain-maintained [nextclade datasets](https://github.com/nextstrain/nextclade_data) for SARS-CoV-2, mpox, influenza A and B, and RSV-A and RSV-B. The organism must be specified as input in the field `organism`, and these align with the nextclade dataset names, i.e. " sars-cov-2", "flu_h1n1pdm_ha", "flu_h1n1pdm_na", "flu_h3n2_ha", "flu_h3n2_na", "flu_vic_ha", "flu_vic_na", "flu_yam_ha", "hMPXV", "hMPXV_B1", "MPXV", "rsv_a" and "rsv_b".

However, nextclade can be used on any organism as long as an an existing, high-quality input reference tree with representative samples on it is provided, in addition to other optional inputs. Contact us if you need help generating your own mutation-annotated tree, or follow the instructions available on the Augur wiki [here](https://docs.nextstrain.org/projects/augur/en/stable/index.html).

!!! info "_Placement_ not _construction_"
    This workflow is not for building a tree from scratch, but rather for the placement of new sequences onto an existing high-quality input reference tree with representative samples on it. In effect, query samples are only compared to reference samples and never to the other query samples.

### Inputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "Sample_to_Ref_Tree"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "Sample_to_Ref_Tree"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///
