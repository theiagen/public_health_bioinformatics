# Workflow Name

## Quick Facts

_Please see the descriptions provided in our documentation contribution guide [here](../contributing/doc_contribution.md#new-page)_

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Link to Workflow Type](../../workflows_overview/workflows_type.md/#link-to-workflow-type) | [Link to Applicable Kingdom](../../workflows_overview/workflows_kingdom.md/#link-to-applicable-kingdom) | PHB <vX.X.X\> | <command-line compatibility\> | <workflow level on terra (set or sample)\> |

## Workflow_Name_On_Terra

_Please provide a description of the workflow._

!!! caption "Workflow Name Diagram"
    _If a workflow diagram is available, add it here. If not, remove this caption. See [the page for TheiaCoV](../workflows/genomic_characterization/theiacov.md#theiacov-workflow-series) for an example._

### Inputs

_Input should be ordered as they appear on Terra. Any additional information regarding the inputs, such as suggestions or guidelines should be found **before** the table. See [the page for Snippy_Streamline](../workflows/phylogenetic_construction/snippy_streamline.md#inputs) for an example._

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| task_name | **variable_name** | Type | Description | Default Value (leave blank if no default) | Required/Optional |

### Workflow Tasks

_Feel free to separate this section into subsections, like "Read QC" and "Alignment" if there are multiple tasks per subsection for easier navigation and readability. See [the page for TheiaMeta](../workflows/genomic_characterization/theiameta.md#workflow-tasks) for an example._

??? task "`tool_name`: Description of tool"
    _Please provide a description of the task and the tool used_

    !!! techdetails "Tool Name Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [link to task on GitHub] |
        | Software Source Code | [link to tool's source code] |
        | Software Documentation | [link to tool's documentation] |
        | Original Publication(s) | [link to tool's publication] |

### Outputs

_Outputs should be in alphabetical order by variable. Any additional information regarding the outputs, such as interpretation suggestions or more details, should be found **after** the table. See [the page for Kraken2](../workflows/standalone/kraken2.md#outputs) for an example._

| **Variable** | **Type** | **Description** |
|---|---|---|
| variable_name | Type | Description |

## References

_Include this section if applicable._

> reference1
<!-- include comments to make sure the `>` blocks are separated -->
> reference2
