# Fetch SRR Accession Workflow

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Data Import](../../workflows_overview/workflows_type.md/#data-import) | [Any Taxa](../../workflows_overview/workflows_kingdom.md/#any-taxa) | PHB v2.3.0 | Yes | Sample-level |

## Fetch SRR Accession

This workflow retrieves the Sequence Read Archive (SRA) accession (SRR) associated with a given sample accession. The primary inputs are BioSample IDs (e.g., SAMN00000000) or SRA Experiment IDs (e.g., SRX000000), which link to sequencing data in the SRA repository.

The workflow uses the fastq-dl tool to fetch metadata from SRA and specifically parses this metadata to extract the associated SRR accession and outputs the SRR accession.

### Inputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="Fetch_SRR_Accession", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Workflow Tasks

This workflow has a single task that performs metadata retrieval for the specified sample accession.

??? task "`fastq-dl`: Fetches SRR metadata for sample accession"
    When provided a BioSample accession or SRA experiment ID, 'fastq-dl' collects metadata and returns the appropriate SRR accession.

    !!! techdetails "fastq-dl Technical Details"
    |  | Links | 
    | --- | --- | 
    | **Task** | [task_fetch_srr_accession.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/data_handling/task_fetch_srr_accession.wdl) |
    | **Software Source Code** | [https://github.com/rpetit3/fastq-dl](https://github.com/rpetit3/fastq-dl) |
    | **Software Documentation** | [https://github.com/rpetit3/fastq-dl#usage](https://github.com/rpetit3/fastq-dl#usage) |

### Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filter_column="Workflow", filter_values="Fetch_SRR_Accession", columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///
