# Find_Shared_Variants

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filter_column="Name", filter_values="[**Find_Shared_Variants**](../workflows/phylogenetic_construction/find_shared_variants.md)", columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level"]) }}

## Find_Shared_Variants_PHB

`Find_Shared_Variants_PHB` is a workflow for concatenating the variant results produced by the `Snippy_Variants_PHB` workflow across multiple samples and reshaping the data to illustrate variants that are shared among multiple samples.

!!! caption "Find_Shared_Variants Workflow Diagram"

    ![Find_Shared_Variants Workflow Diagram](../../assets/figures/Find_Shared_Variants_PHB.png)

### Inputs

The primary intended input of the workflow is the `snippy_variants_results` output from `Snippy_Variants_PHB` or the `theiaeuk_snippy_variants_results` output of the TheiaEuk workflow. Variant results files from other tools may not be compatible at this time.

All variant data included in the sample set should be generated from aligning sequencing reads to the **same reference genome**. If variant data was generated using different reference genomes, shared variants cannot be identified and results will be less useful.

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "Find_Shared_Variants"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Tasks

{{ include_md("common_text/shared_variants_task.md", condition="find_shared_variants") }}

### Outputs

The outputs of this workflow are the `concatenated_variants` file and the `shared_variants_table` file.

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "Find_Shared_Variants"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///
