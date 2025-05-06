# Rename_FASTQ

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Standalone](../../workflows_overview/workflows_type.md/#standalone) | [Any Taxa](../../workflows_overview/workflows_kingdom.md/#any-taxa) | PHB v2.1.0 | Yes | Sample-level |

## Rename_FASTQ_PHB

This sample-level workflow receives a read file or a pair of read files (FASTQ), compressed or uncompressed, and returns a new, renamed and compressed FASTQ file.

### Inputs

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="Rename_FASTQ", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Outputs

If a reverse read (`read2`) is provided, the files get renamed to the provided `new_filename` input with the notation `<new_filename>_R1.fastq.gz` and `<new_filename>_R2.fastq.gz`. If only `read1` is provided, the file is renamed to `<new_filename>.fastq.gz`.

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/tables/all_outputs.tsv", input_table=False, filter_column="Workflow", filter_values="Rename_FASTQ", columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///
