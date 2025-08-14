# NCBI_Scrub

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Name": "[**NCBI_Scrub**](../workflows/standalone/ncbi_scrub.md)"}, columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level", "Dockstore"]) }}

## NCBI Scrub Workflows

NCBI Scrub, also known as the human read removal tool (HRRT), is based on the [SRA Taxonomy Analysis Tool](https://doi.org/10.1186/s13059-021-02490-0) that will take as input a FASTQ file, and produce as output a FASTQ file in which all reads identified as potentially of human origin are either removed (default) or masked with 'N'.
There are three Kraken2 workflows:

- `NCBI_Scrub_PE` is compatible with **Illumina paired-end data**
- `NCBI_Scrub_SE` is compatible with **Illumina single-end data**

### Inputs

!!! caption ""
    === "NCBI_Scrub_PE"
        /// html | div[class="searchable-table"]

        {{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "NCBI_Scrub_PE"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"], indent=8) }}
        ///

    === "NCBI_Scrub_SE"
        /// html | div[class="searchable-table"]

        {{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "NCBI_Scrub_SE"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"], indent=8) }}
        ///

### Workflow Tasks

This workflow is composed of two tasks, one to dehost the input reads and another to screen the clean reads with kraken2 and the viral+human database.

{{ include_md("common_text/ncbi_scrub_task.md") }}
{{ include_md("common_text/kraken2_task.md", condition="theiacov") }}

### Outputs

!!! caption ""
    === "NCBI_Scrub_PE"
        /// html | div[class="searchable-table"]

        {{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "NCBI_Scrub_PE"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"], indent=8) }}

        ///

    === "NCBI_Scrub_SE"
        /// html | div[class="searchable-table"]

        {{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "NCBI_Scrub_SE"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"], indent=8) }}
        
        ///
