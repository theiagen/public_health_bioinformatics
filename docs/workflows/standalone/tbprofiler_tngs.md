# TBProfiler_tNGS

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Standalone](../../workflows_overview/workflows_type.md/#standalone) | [Bacteria](../../workflows_overview/workflows_kingdom.md/#bacteria) | PHB v3.0.1 | Yes | Sample-level |

## TBProfiler_tNGS_PHB

This workflow is still in experimental research stages. Documentation is minimal as changes may occur in the code; it will be fleshed out when a stable state has been achieved.

### Inputs

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/input_tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="TBProfiler_tNGS", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Terra Outputs

<div class="searchable-table" markdown="1">

| **Variable** | **Type** | **Description** |
|---|---|---|
| tbp_parser_average_genome_depth | Float | The mean depth of coverage across all target regions included in the analysis |
| tbp_parser_coverage_report | File | A file containing the breadth of coverage across each target loci |
| tbp_parser_docker | String | The docker image and version tag for the tbp_parser tool |
| tbp_parser_genome_percent_coverage | Float | The percent breadth of coverage across the entire genome  |
| tbp_parser_laboratorian_report_csv | File | An output file containing information regarding each mutation and its associated drug resistance profile in a CSV file. This file also contains two interpretation fields -- "Looker" and "MDL" which are generated using the CDC's expert rules for interpreting the severity of potential drug resistance mutations. |
| tbp_parser_lims_report_csv | File | An output file formatted specifically for STAR LIMS. This CSV report summarizes the highest severity mutations for each antimicrobial and lists the relevant mutations for each gene. |
| tbp_parser_looker_report_csv | File | An output file that contains condensed information suitable for generating a dashboard in Google's Looker studio. |
| tbp_parser_version | String | The version number of tbp_parser |
| tbprofiler_dr_type | String | The drug resistance category as determined by TBProfiler |
| tbprofiler_main_lineage | String | The Mycobacterium tuberculosis lineage assignment as made by TBProfiler |
| tbprofiler_median_coverage | Int | The median depth of coverage across the target loci |
| tbprofiler_num_dr_variants | String | The total number of drug resistance conferring variants detected by TBProfiler |
| tbprofiler_num_other_variants | String | The total number of non-drug resistance conferring variants detected by TBProfiler |
| tbprofiler_output_alignment_bai | File | The index file associated with the binary alignment map of the input reads against the H37Rv genome |
| tbprofiler_output_alignment_bam | File | The binary alignment map of the input reads against the H37Rv genome |
| tbprofiler_pct_reads_mapped | Float | The percentage of reads that successfully mapped to the H37Rv genome |
| tbprofiler_report_csv | File | The raw output file from TBProfiler |
| tbprofiler_report_json | File | The json output file from TBProfiler |
| tbprofiler_report_tsv | File | The TSV output file from TBProfiler |
| tbprofiler_resistance_genes | String | The genes in which a mutation was detected that may be resistance conferring |
| tbprofiler_sub_lineage | String | The Mycobacterium tuberculosis sub-lineage assignment as made by TBProfiler |
| tbprofiler_tngs_wf_analysis_date | String | The date on which the workflow was run |
| tbprofiler_tngs_wf_version | String | The version of the tbprofiler_tngs workflow used for this analysis |
| tbprofiler_version | String | The version of TBProfiler used for this analysis |
| trimmomatic_docker | String | The docker image used for the trimmomatic module in this workflow |
| trimmomatic_read1_trimmed | File | The read1 file post trimming |
| trimmomatic_read2_trimmed | File | The read2 file post trimming |
| trimmomatic_stats | File | The read trimming statistics |
| trimmomatic_version | String | The version of trimmomatic used in this analysis |

</div>
