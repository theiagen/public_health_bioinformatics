# Concatenate Illumina Lanes

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Standalone](../../workflows_overview/workflows_type.md/#standalone) | [Any Taxa](../../workflows_overview/workflows_kingdom.md/#any-taxa) | PHB 2.3.0 | Yes | Sample-level |

## Concatenate_Illumina_Lanes_PHB

Some Illumina machines produce multi-lane FASTQ files for a single sample. This workflow concatenates the multiple lanes into a single FASTQ file per read type (forward or reverse).

### Inputs

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/input_tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="Concatenate_Illumina_Lanes", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Workflow Tasks

This workflow concatenates the Illumina lanes for forward and (if provided) reverse reads. The output files are named as followed:

- Forward reads: `<samplename>_merged_R1.fastq.gz`
- Reverse reads: `<samplename>_merged_R2.fastq.gz`

### Outputs

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/tables/all_outputs.tsv", input_table=False, filter_column="Workflow", filter_values="Concatenate_Illumina_Lanes", columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///
