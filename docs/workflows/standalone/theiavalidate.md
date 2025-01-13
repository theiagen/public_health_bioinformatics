# TheiaValidate

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Standalone](../../workflows_overview/workflows_type.md/#standalone) | [Any Taxa](../../workflows_overview/workflows_kingdom.md/#any-taxa) | PHB v2.0.0 | No | |

## TheiaValidate_PHB

!!! caption "TheiaValidate Workflow Diagram"
    ![TheiaValidate Workflow Diagram](../../assets/figures/TheiaValidate.png)

TheiaValidate performs basic comparisons between user-designated columns in two separate tables. We anticipate this workflow being run to determine if any differences exist between version releases or two workflows, such as TheiaProk_ONT vs TheiaProk_Illumina_PE. A summary PDF report is produced in addition to a Excel spreadsheet that lists the values for any columns that do not have matching content for a sample.

!!! warning
    The two tables being compared **must** have both identical sample names and an equal number of samples. If not, validation will not work or (in the case of unequal number of samples) not be attempted.

In order to enable this workflow to function for different workflow series, we require users to provide a list of columns they want to compare between the two tables. Feel free to use the information below that Theiagen uses to compare versions of the three main workflow series as a _**starting point**_ for your own validations:

!!! tool "Validation Starting Points"
    | Workflow Series | Validation Criteria TSV | Columns to Compare |
    |---|---|---|
    | TheiaCoV Workflows | [TheiaCov Validation Criteria](../../assets/files/theiavalidate/theiacov-validation-criteria.txt) | abricate_flu_subtype,abricate_flu_type,assembly_length_unambiguous,assembly_mean_coverage,irma_subtype,irma_type,kraken_human,kraken_human_dehosted,kraken_sc2,kraken_sc2_dehosted,kraken_target_org,kraken_target_org_dehosted,nextclade_aa_dels,nextclade_aa_subs,nextclade_clade,nextclade_lineage,nextclade_tamiflu_resistance_aa_subs,num_reads_clean1,num_reads_clean2,number_N,pango_lineage,percent_reference_coverage,vadr_num_alerts |
    | TheiaEuk Workflows | [TheiaEuk Validation Criteria](../../assets/files/theiavalidate/theiaeuk-validation-criteria.txt) | assembly_length,busco_results,clade_type,est_coverage_clean,est_coverage_raw,gambit_predicted_taxon,n50_value,num_reads_clean1,num_reads_clean2,number_contigs,quast_gc_percent,theiaeuk_snippy_variants_hits |
    | TheiaProk Workflows | [TheiaProk Validation Criteria](../../assets/files/theiavalidate/theiaprok-validation-criteria.txt) | abricate_abaum_plasmid_type_genes,agrvate_agr_group,amrfinderplus_amr_core_genes,amrfinderplus_amr_plus_genes,amrfinderplus_stress_genes,amrfinderplus_virulence_genes,ani_highest_percent,ani_top_species_match,assembly_length,busco_results,ectyper_predicted_serotype,emmtypingtool_emm_type,est_coverage_clean,est_coverage_raw,gambit_predicted_taxon,genotyphi_final_genotype,hicap_genes,hicap_serotype,kaptive_k_type,kleborate_genomic_resistance_mutations,kleborate_key_resistance_genes,kleborate_mlst_sequence_type,legsta_predicted_sbt,lissero_serotype,meningotype_serogroup,midas_primary_genus,midas_secondary_genus,midas_secondary_genus_abundance,n50_value,ngmaster_ngmast_sequence_type,ngmaster_ngstar_sequence_type,num_reads_clean1,num_reads_clean2,number_contigs,pasty_serogroup,pbptyper_predicted_1A_2B_2X,plasmidfinder_plasmids,poppunk_gps_cluster,seqsero2_predicted_serotype,seroba_ariba_serotype,seroba_serotype,serotypefinder_serotype,shigatyper_ipaB_presence_absence,shigatyper_predicted_serotype,shigeifinder_cluster,shigeifinder_serotype,sistr_predicted_serotype,sonneityping_final_genotype,spatyper_type,srst2_vibrio_serogroup,staphopiasccmec_types_and_mecA_presence,tbprofiler_main_lineage,tbprofiler_resistance_genes,ts_mlst_predicted_st,virulencefinder_hits |

If additional validation metrics are desired, the user has the ability to provide a `validation_criteria_tsv` file that specifies what type of comparison should be performed. There are several options for additional validation checks:

- **EXACT** performs an exact string match and counts the number of exact match failures/differences
- **IGNORE** does not check the values and says there are 0 failures
- **SET** checks list items (such as `amrfinder_plus_genes` which is a comma-delimited list of genes) for identical content — order does not matter; that is, `mdsA,mdsB` is determined to be same as `mdsB,mdsA`. The EXACT match does not consider these to be the same, but the SET match does.
- **<PERCENT_DIFF\>**, which is an actual decimal value such as **0.02**, calculates the percent difference between _numerical_ columns. If the columns are not numerical, this function will **not** work and will lead to workflow failure. For example, if the decimal percentage is 0.02, the test will indicate a failure if the values in the two columns are more than 2% different.
- Dates, integers, and object-type values are ignored and indicate 0 failures.

### File Comparisons

If a column consists of only GCP URIs (Google Cloud file paths), the files will be localized and compared with either an EXACT match or a SET match. In the SET match, the lines in the file are ordered before comparison. Results are returned to the summary table as expected. The results of each file comparison can be found in the `theiavalidate_diffs`  output column.

### Inputs

<div class="searchable-table" markdown="1">

Please note that all string inputs **must** be enclosed in quotation marks; for example, "column1,column2" or "workspace1"

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| theiavalidate | **columns_to_compare** | String | A comma-separated list of the columns the user wants to compare. Do not include whitespace. |  | Required |
| theiavalidate | **output_prefix** | String | The prefix for the output files |  | Required |
| theiavalidate | **table1_name** | String | The name of the first table |  | Required |
| theiavalidate | **table2_name** | String | The name of the second table |  | Required |
| theiavalidate | **terra_project1_name** | String | The name of the Terra project where table1_name can be found |  | Required |
| theiavalidate | **terra_workspace1_name** | String | The name of the Terra workspace where table1_name can be found |  | Required |
| theiavalidate | **column_translation_tsv** | File | If the user wants to link two columns of different names, they may supply a TSV file that provides a "column translation" between the two files (see the section below this table). |  | Optional |
| theiavalidate | **terra_project2_name** | String | If the table2_name is located in a different Terra project, indicate it here. Otherwise, the workflow will look for table2_name in the Terra project indicated in terra_project1_name. | value for `terra_project1_name` | Optional |
| theiavalidate | **terra_workspace2_name** | String | If the table2_name is located in a different Terra workspace, indicate it here. Otherwise, the workflow will look for table2_name in the Terra workspace indicated in terra_workspace1_name. | value for `terra_workspace1_name` | Optional |
| theiavalidate | **validation_criteria_tsv** | File | If the user wants to specify a different comparison than the default exact string match, they may supply a TSV file that indicates the different options (see the section below this table). |  | Optional |
| compare_two_tsvs | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| compare_two_tsvs | **debug_output** | Boolean | Set to true to enable more outputs; useful when debugging | FALSE | Optional |
| compare_two_tsvs | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| compare_two_tsvs | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/theiavalidate:0.1.0 | Optional |
| compare_two_tsvs | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 4 | Optional |
| compare_two_tsvs | **na_values** | String | If the user knows a particular value in either table that they would like to be considered N/A, they can indicate those values in a comma-separated list here. Any changes here will overwrite the default and not append to the default list. Do not include whitespace. | -1.#IND,1.#QNAN,1.#IND,-1.#QNAN,#N/A,N/A,n/a,,#NA,NULL,null,NaN,-NaN,nan,-nan,None | Optional |
| export_two_tsvs | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional |
| export_two_tsvs | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 10 | Optional |

</div>

The optional `validation_criteria_tsv` file takes the following format (tab-delimited; _a header line is required_):

```text linenums="1"
column_name	criteria
columnB SET
columnC	IGNORE
columnD	0.01
columnE	EXACT
```

Please see above for a description of all available criteria options (EXACT, IGNORE, SET, <PERCENT_DIFF>).

The optional `column_translation_tsv` file takes the following format (tab-delimited; _there can be **no** header line_):

```text linenums="1"
column_name_in_table1	column_name_in_table2
column_name_in_table2	column_name_in_table1
internal_column_name	display_column_name
```

Please note that the name in the **second column** will be displayed and used in all output files.

!!! warning "Known Bug"
    There must be _**more**_ than one line in the `column_translation_tsv` file or else this error will appear: `AttributeError: 'str' object has no attribute 'to_dict'`. To fix this error, add an additional line in the `column_translation_tsv` file, like the following: `columnA	columnA`

!!! warning "Known Bug"
    If performing a <PERCENT_DIFF> comparison, all samples must have values for that column.

!!! info "Call Caching Disabled"
    If using TheiaValidate workflow version 1.3.0 or higher, the call-caching feature of Terra has been DISABLED to ensure that the workflow is run from the beginning and data is compared fresh. Call-caching will not be enabled, even if the user checks the box ✅ in the Terra workflow interface.

### Outputs

<div class="searchable-table" markdown="1">

| **Variable** | **Type** | **Description** |
|---|---|---|
| theiavalidate_criteria_differences | File | A TSV file that lists only the differences that fail to meet the validation criteria |
| theiavalidate_date | String | The date the analysis was run |
| theiavalidate_diffs | Array[File] | An array of files with a single file for each file comparison performed; only has values if a column with files is compared |
| theiavalidate_exact_differences | File | A TSV file that lists all exact string match differences between samples |
| theiavalidate_filtered_input_table1 | File | The first data table used for validation after removing unexamined columns and translating column names |
| theiavalidate_filtered_input_table2 | File | The second data table used for validation after removing unexamined columns and translating column names |
| theiavalidate_report | File | A PDF summary report |
| theiavalidate_status | String | Indicates whether or not validation was attempted |
| theiavalidate_version | String | The version of the TheiaValidate Python Docker |
| theiavalidate_wf_version | String | The version of the PHB repository |

</div>

### Example Data and Outputs

To help demonstrate how TheiaValidate works, please observe the following example and outputs:

???+ toggle "Table1"
    | entity:example_table1_id | columnA-string | columnB-set | columnC-ignore | columnD-float | columnE-missing |
    | --- | --- | --- | --- | --- | --- |
    | sample1 | option1 | item1,item2,item3 | cheese | 1000 | present |
    | sample2 | option1 | item1,item3,item2 | cheesecake | 12 | present |
    | sample3 | option2 | item1,item2,item3 | cake | 14 | present |
    | sample4 | option1 | item2,item1 | cakebatter | 3492 |  |
    | sample5 | option2 | item1,item2 | batter | 3 | present |
  
???+ toggle "Table2"
    | entity:example_table2_id | columnA-string | columnB-set | columnC-ignore | columnD-float | missing |
    | --- | --- | --- | --- | --- | --- |
    | sample1 | option1 | item1,item3,item2 | cheesecake | 999 | present |
    | sample2 | option2 | item1,item2,item3 | batter | 12 | present |
    | sample3 | option1 | item1,item2 | cheese | 24 |  |
    | sample4 | option1 | item1,item2 | cakebatter | 728 |  |
    | sample5 | option2 | item1,item2,item3 | batter | 4 | present |

???+ toggle "Validation Criteria"
    | column | criteria |
    | --- | --- |
    | columnB-set | SET |
    | columnC-ignore | IGNORE |
    | columnD-float | 0.01 |
    | columnE-missing | EXACT |

???+ toggle "Column Translation"
    | missing | columnE-missing |
    | --- | --- |
    | columnA-string | columnA-string |

    _Note: the second row translating_ `columnA-string` _to itself is included to prevent the known bug explained above._

If the above inputs are provided, then the following output files will be generated:

[filtered_example_table1.tsv](../../assets/files/theiavalidate/filtered_example_table1.tsv)

[filtered_example_table2.tsv](../../assets/files/theiavalidate/filtered_example_table2.tsv)

[example_summary.pdf](../../assets/files/theiavalidate/example_summary.pdf)

[example_exact_differences.tsv](../../assets/files/theiavalidate/example_exact_differences.tsv)

[example_validation_criteria_differences.tsv](../../assets/files/theiavalidate/example_validation_criteria_differences.tsv)