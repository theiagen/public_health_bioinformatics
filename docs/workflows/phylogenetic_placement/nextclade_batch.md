# Nextclade_Batch

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Name": "[**Nextclade_Batch**](../workflows/phylogenetic_placement/nextclade_batch.md)"}, columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level"]) }}

## Nextclade_Batch_PHB

Nextclade Batch rapidly calls mutations, places samples on a reference phylogenetic tree, and rapidly genotypes batches of samples using [Nextclade](https://docs.nextstrain.org/projects/nextclade/en/stable/index.html). Phylogenetic placement is done by comparing the mutations of the query sequence (relative to the reference) with the mutations of every node and tip in the reference tree, and finding the node which has the most similar set of mutations. This operation is repeated for each query sequence, until all of them are placed onto the tree. This workflow uses the Nextstrain-maintained [nextclade datasets](https://github.com/nextstrain/nextclade_data) for manually inputted datasets or downloaded datasets (e.g. SARS-CoV-2, mpox, influenza A and B, HIV, and RSV-A and RSV-B).

Contact us if you need help generating your own mutation-annotated tree, or follow the instructions available on the Augur wiki [here](https://docs.nextstrain.org/projects/augur/en/stable/index.html).

!!! info "_Placement_ not _construction_"
    This workflow is not for building a tree from scratch, but rather for genotyping and placement of new sequences onto an existing high-quality input reference tree with representative samples on it. In effect, query samples are only compared to reference samples and never to the other query samples.

### Inputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "Nextclade_Batch"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "Nextclade_Batch"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///
