# Workflow Name

## Quick Facts

_Use the following render_tsv_table macro call that is provided in code below to generate the **quick facts table**, replacing the fields marked with `<>` with the appropriate values. Please note that the macro_ **result** _is seen on the web browser, not the macro call itself._

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Name": "[**<your-workflow-name\>**](../workflows/<your-workflow-type>/<your-workflow-name>.md)"}, columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level", "Dockstore"]) }}

## Workflow_Name_On_Terra

_Please provide a description of the workflow in complete sentences._

!!! caption "Workflow Name Diagram"
    _If a workflow diagram is available, add it here. If not, remove this caption admonition. See [the page for TheiaCoV](../workflows/genomic_characterization/theiacov.md#theiacov-workflow-series) for an example._

### Inputs

_Any additional information regarding the inputs, such as suggestions or guidelines should be found **before** the table. See [the page for Snippy_Streamline](../workflows/phylogenetic_construction/snippy_streamline.md#inputs) for an example._

_Add your inputs anywhere in the `assets/tables/all_inputs.tsv` table. If your workflow has the same type of input as a different workflow, check the inputs table to see if that input already exists. If so, add your workflow name to the comma-delimited list for that input in the appropriate column._

_Use the following macro call to generate the inputs table. Adjust the fields marked with `<>` with the appropriate values.  Please note that the macro_ **result** _is seen on the web browser, not the macro call itself._

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "<your-workflow-name\>"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Workflow Tasks

_Feel free to separate this section into subsections, like "Read QC" and "Alignment" if there are multiple tasks per subsection for easier navigation and readability. See [the page for TheiaMeta](../workflows/genomic_characterization/theiameta.md#workflow-tasks) for an example._

_If your workflow uses a task that is modular and can be used in other contexts, please add that information to the `docs/common_text` directory in a new page (see a template in the `/common_text/template_task.md` file) and use the following macro call to include it here. Adjust the fields marked with `<>` with the appropriate values. Please note that the macro_ **result** _is seen on the web browser, not the macro call itself._

{{ include_md("common_text/template_task.md", condition="condition") }}

### Outputs

_Add your outputs to the `assets/tables/all_outputs.tsv` table. Location doesn't matter. If your workflow has the same type of output as a different workflow, check the output table to see if that output already exists. If so, add your workflow name to the comma-delimited list in the appropriate column._

_Use the following macro call to generate the outputs table. Adjust the fields marked with `<>` with the appropriate values. Please note that the macro_ **result** _is seen on the web browser, not the macro call itself._

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "<your-workflow-name\>"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///

_Any additional information regarding the outputs, such as interpretation suggestions or more details, should be found **after** the table. See [the page for Kraken2](../workflows/standalone/kraken2.md#outputs) for an example._

## References

_Include this section if applicable._

> reference1
<!-- include comments to make sure the `>` blocks are separated -->
> reference2
