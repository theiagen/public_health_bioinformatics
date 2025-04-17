# TBProfiler_tNGS

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Standalone](../../workflows_overview/workflows_type.md/#standalone) | [Bacteria](../../workflows_overview/workflows_kingdom.md/#bacteria) | PHB vX.X.X | Yes | Sample-level |

## TBProfiler_tNGS_PHB

This workflow is still in experimental research stages. Documentation is minimal as changes may occur in the code; it will be fleshed out when a stable state has been achieved.

### Inputs

<div class="searchable-table" markdown="1">

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| tbprofiler_tngs | **read1** | File | Illumina forward read file in FASTQ file format (compression optional) |  | Required |
| tbprofiler_tngs | **read2** | File | Illumina reverse read file in FASTQ file format (compression optional) |  | Required |
| tbprofiler_tngs | **samplename** | String | Name of sample to be analyzed |  | Required |
| tbp_parser | **config** | File | The configuration file to use, in YAML format (overrides all other arguments except input_json and input_bam) |  | Optional |
| tbp_parser | **coverage_regions_bed** | File | A file that contains the regions to perform coverage analysis on |  | Optional |
| tbp_parser | **min_percent_coverage** | Float | The minimum percentage of a region to exceed the minimum depth for a region to pass QC in tbp_parser | 100.0 | Optional |
| tbp_parser | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional |
| tbp_parser | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| tbp_parser | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/tbp-parser:2.4.5 | Optional |
| tbp_parser | **etha237_frequency** | Float | Minimum frequency for a mutation in ethA at protein position 237 to pass QC in tbp-parser | 0.1 | Optional |
| tbp_parser | **expert_rule_regions_bed** | File | A file that contains the regions where R mutations and expert rules are applied |  | Optional |
| tbp_parser | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 4 | Optional |
| tbp_parser | **min_depth** | Int | Minimum depth for a variant to pass QC in tbp_parser | 10 | Optional |
| tbp_parser | **min_frequency** | Float | Minimum allele frequency for a variant to pass QC in tbp-parser | 0.1 | Optional |
| tbp_parser | **min_read_support** | Int | Minimum read support for a variant to pass QC in tbp-parser | 10 | Optional |
| tbp_parser | **operator** | String | Fills the "operator" field in the tbp_parser output files |  | Optional |
| tbp_parser | **rpob449_frequency** | Float | Minimum frequency for a mutation at protein position 449 to pass QC in tbp-parser | 0.1 | Optional |
| tbp_parser | **rrl_frequency** | Float | Minimum frequency for a mutation in rrl to pass QC in tbp-parser | 0.1 | Optional |
| tbp_parser | **rrl_read_support** | Int | Minimum read support for a mutation in rrl to pass QC in tbp-parser | 10 | Optional |
| tbp_parser | **rrs_frequency** | Float | Minimum frequency for a mutation in rrs to pass QC in tbp-parser | 0.1 | Optional |
| tbp_parser | **rrs_read_support** | Int | Minimum read support for a mutation in rrs to pass QC in tbp-parser | 10 | Optional |
| tbp_parser | **sequencing_method** | String | Fills out the "seq_method" field in the tbp_parser output files |  | Optional |
| tbp_parser | **tbp_parser_debug** | Boolean | Activate the debug mode on tbp_parser; increases logging outputs | FALSE | Optional |
| tbprofiler | **cov_frac_threshold** | Int | A cutoff used to calculate the fraction of the region covered by ≤ this value | 1 | Optional |
| tbprofiler | **cpu** | Int | Number of CPUs to allocate to the task | 8 | Optional |
| tbprofiler | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| tbprofiler | **mapper** | String | The mapping tool used in TBProfiler to align the reads to the reference genome; see TBProfiler's original documentation for available options. | bwa | Optional |
| tbprofiler | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 16 | Optional |
| tbprofiler | **min_af** | Float | The minimum allele frequency to call a variant | 0.1 | Optional |
| tbprofiler | **min_af_pred** | Float | The minimum allele frequency to use a variant for resistance prediction | 0.1 | Optional |
| tbprofiler | **min_depth** | Int | The minimum depth for a variant to be called. | 10 | Optional |
| tbprofiler | **ont_data** | Boolean | Internal component; do not modify |  | Do not modify, Optional |
| tbprofiler | **tbprofiler_custom_db** | File | TBProfiler uses by default the TBDB database; if you have a custom database you wish to use, you must provide a custom database in this field and set tbprofiler_run_custom_db to true |  | Optional |
| tbprofiler | **tbprofiler_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/tbprofiler:6.6.3 | Optional |
| tbprofiler | **tbprofiler_run_custom_db** | Boolean |  | FALSE | Optional |
| tbprofiler | **variant_caller** | String | Select a different variant caller for TBProfiler to use by writing it in this block; see TBProfiler's original documentation for available options. | freebayes | Optional |
| tbprofiler | **variant_calling_params** | String | Enter additional variant calling parameters in this free text input to customize how the variant caller works in TBProfiler |  | Optional |
| tbprofiler | **bases_to_crop** | Int | Indicate the number of bases to remove from the start and end of the read | 0 | Optional |
| trimmomatic_pe | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional |
| trimmomatic_pe | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| trimmomatic_pe | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/trimmomatic:0.39 | Optional |
| trimmomatic_pe | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| trimmomatic_pe | **trimmomatic_args** | String | Additional arguments to pass to trimmomatic. "-phred33" specifies the Phred Q score encoding which is almost always phred33 with modern sequence data. | -phred33 | Optional |
| trimmomatic_pe | **trimmomatic_min_length** | Int | Specifies minimum length of each read after trimming to be kept | 75 | Optional |
| trimmomatic_pe | **trimmomatic_quality_trim_score** | Int | The trimming quality score | 30 | Optional |
| trimmomatic_pe | **trimmomatic_window_size** | Int | The window size for trimming | 4 | Optional |
| version_capture | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional |

</div>

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
