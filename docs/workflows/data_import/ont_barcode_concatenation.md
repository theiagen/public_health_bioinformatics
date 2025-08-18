# ONT_Barcode_Concatenation

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Name": "[**ONT_Barcode_Concatenation**](../workflows/data_import/ont_barcode_concatenation.md)"}, columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level", "Dockstore"]) }}

## ONT_Barcode_Concatenation_PHB

Unconcatenated ONT data is the bane of all humanity. This workflow will automatically concatenate all reads in a given folder and upload those reads to a Terra data table.

We recommend running this workflow with **"Run workflow with inputs defined by file paths"** selected in Terra. This will allow you to upload your data files and provide the necessary information for the workflow without having to specify a data table. There are no outputs for this workflow, as the data is added to either a new or existing table in your workspace.

### Inputs

!!! warning "Default Behavior"

NESTED STUFF

#### Uploading unconcatenated ONT reads to Terra {% raw %} {#data-upload} {% endraw %}

Using the Terra data uploader is **not recommended**.

#### Finding the `input_bucket_path` and `output_bucket_path` {% raw %} {#file-paths} {% endraw %}

You can find the file paths by ....

#### Creating a `barcode_renaming_file` {% raw %} {#barcode-renaming} {% endraw %}

Make a file!

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "ONT_Barcode_Concatenation"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Workflow Tasks

_Feel free to separate this section into subsections, like "Read QC" and "Alignment" if there are multiple tasks per subsection for easier navigation and readability. See [the page for TheiaMeta](../workflows/genomic_characterization/theiameta.md#workflow-tasks) for an example._

_If your workflow uses a task that is modular and can be used in other contexts, please add that information to the `docs/common_text` directory in a new page (see a template in the `/common_text/template_task.md` file) and use the following macro call to include it here. Adjust the fields marked with `<>` with the appropriate values. Please note that the macro_ **result** _is seen on the web browser, not the macro call itself._

{{ include_md("common_text/template_task.md", condition="condition") }}

### Outputs

Your concatenated ONT data will automatically appear in your workspace in the table of choice with information in the following four fields:

- Sample name (under the `terra_table_name`_id column), which will be either the name of the parent folder or the remapped name indicated by the `barcode_renaming_file` input.
- The concatenated ONT data in the `read1` column
- The name of the workflow (`ONT_Barcode_Concatenation_PHB`) under the `table_created_by` column, to indicate the samples were added by this workflow.
- The date of upload/when the workflow was run under the `upload_date` column
