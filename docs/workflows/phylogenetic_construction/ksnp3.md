# kSNP3

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Phylogenetic Construction](../../workflows_overview/workflows_type.md/#phylogenetic-construction) | [Bacteria](../../workflows_overview/workflows_kingdom.md/#bacteria), [Mycotics](../../workflows_overview/workflows_kingdom.md#mycotics), [Viral](../../workflows_overview/workflows_kingdom.md/#viral) | PHB v3.0.0 | Yes; some optional features incompatible | Set-level |

## kSNP3_PHB

The kSNP3 workflow is for phylogenetic analysis of bacterial genomes using single nucleotide polymorphisms (SNPs). The kSNP3 workflow identifies SNPs amongst a set of genome assemblies, then calculates a number of phylogenetic trees based on those SNPs:

- **Pan-genome phylogenetic trees:** The term "pan-genome" is used here to describe the collective genetic content amongst the set of genomes, including regions outside of genes and other coding sequences.  Outputs based on the pan-genome are labeled with `_pan`.
- **Core-genome phylogenetic trees:** The kSNP3 workflow will also generate phylogenetic trees based on the core genome (genetic content that is present in all members of the set of genomes). Outputs based on the core-genome are labeled with `_core`.

This workflow also features an optional module, `summarize_data` that creates a presence/absence matrix for the analyzed samples from a list of indicated columns (such as AMR genes, plasmid types etc.). If the `phandango_coloring` variable is set to `true`, this will be formatted for visualization in [Phandango](https://jameshadfield.github.io/phandango/#/), else it can be viewed in Excel.

You can learn more about the kSNP3 workflow, including how to visualize the outputs with MicrobeTrace in the following video: **📺 [Using KSNP3 in Terra and Visualizing Bacterial Genomic Networks in MicrobeTrace](https://www.youtube.com/watch?v=iRpNDun46R8)**

### Inputs

<div class="searchable-table" markdown="1">

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| ksnp3_workflow | **assembly_fasta** | Array[File] | The assembly files to be analyzed | | Required |
| ksnp3_workflow | **cluster_name** | String | Free text string used to label output files | | Required |
| ksnp3_workflow | **samplename** | Array[String] | The set of sample names | | Required |
| core_ksnp3_shared_snps_task | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| core_reorder_matrix | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| core_reorder_matrix | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| core_reorder_matrix | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/mykrobe:0.12.1 | Optional |
| core_reorder_matrix | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional |
| core_snp_dists | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional |
| core_snp_dists | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 50 | Optional |
| core_snp_dists | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/snp-dists:0.8.2 | Optional |
| core_snp_dists | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional |
| ksnp3_task | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional |
| ksnp3_task | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| ksnp3_task | **docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/ksnp3:3.1 | Optional |
| ksnp3_task | **kmer_size** | Int | The length of kmer containing the SNP you want kSNP3 to use | 19 | Optional |
| ksnp3_task | **ksnp3_args** | String | Additional arguments you want kSNP3 to use; e.g., "-ML" or "-NJ" |  | Optional |
| ksnp3_task | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| ksnp3_task | **previous_ksnp3_snps** | File | File with existing SNPs for the current run to be appended to.  |  | Optional |
| ksnp3_workflow | **data_summary_column_names** | String | A comma-separated list of the column names from the sample-level data table for generating a data summary (presence/absence .csv matrix); e.g., "amrfinderplus_amr_genes,amrfinderplus_virulence_genes" |  | Optional |
| ksnp3_workflow | **data_summary_terra_project** | String | The billing project for your current workspace. This can be found after the "#workspaces/" section in the workspace's URL |  | Optional |
| ksnp3_workflow | **data_summary_terra_table** | String | The name of the sample-level Terra data table that will be used for generating a data summary |  | Optional |
| ksnp3_workflow | **data_summary_terra_workspace** | String | The name of the Terra workspace you are in. This can be found at the top of the webpage, or in the URL after the billing project. |  | Optional |
| ksnp3_workflow | **midpoint_root_tree** | Boolean | If true, midpoint root the final tree | FALSE | Optional |
| ksnp3_workflow | **phandango_coloring** | Boolean | Boolean variable that tells the data summary task and the reorder matrix task to include a suffix that enables consistent coloring on Phandango; by default, this suffix is not added. To add this suffix set this variable to true. | FALSE | Optional |
| pan_reorder_matrix | **cpu** | Int | Number of CPUs to allocate to the task | 100 | Optional |
| pan_reorder_matrix | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 2 | Optional |
| pan_reorder_matrix | **docker** | String | The Docker container to use for the task | 100 | Optional |
| pan_reorder_matrix | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | us-docker.pkg.dev/general-theiagen/staphb/mykrobe:0.12.1 | Optional |
| pan_snp_dists | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional |
| pan_snp_dists | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 50 | Optional |
| pan_snp_dists | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/snp-dists:0.8.2 | Optional |
| pan_snp_dists | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional |
| summarize_data | **cpu** | Int | Number of CPUs to allocate to the task | 8 | Optional |
| summarize_data | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| summarize_data | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16 | Optional |
| summarize_data | **id_column_name** | String | If the sample IDs are in a different column to samplenames, it can be passed here and it will be used instead. |  | Optional |
| summarize_data | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| version_capture | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional |

</div>

### Workflow Actions

The `ksnp3` workflow is run on the set of assembly files to produce both pan-genome and core-genome phylogenies. This also results in alignment files which - are used by [`snp-dists`](https://github.com/tseemann/snp-dists) to produce a pairwise SNP distance matrix for both the pan-genome and core-genomes.

If you fill out the `data_summary_*` and `sample_names` optional variables, you can use the optional `summarize_data` task. The task takes a comma-separated list of column names from the Terra data table, which should each contain a list of comma-separated items. For example, `"amrfinderplus_virulence_genes,amrfinderplus_stress_genes"` (with quotes, comma separated, no spaces) for these output columns from running TheiaProk. The task checks whether those comma-separated items are present in each row of the data table (sample), then creates a CSV file of these results. The CSV file indicates presence (TRUE) or absence (empty) for each item. By default, the task adds a Phandango coloring tag to group items from the same column, but you can turn this off by setting `phandango_coloring` to `false`.

??? toggle "**Example output CSV**"

    ```text linenums="1"
    Sample_Name,aph(3')-IIa,blaCTX-M-65,blaOXA-193,tet(O)
    sample1,TRUE,,TRUE,TRUE
    sample2,,,FALSE,TRUE
    sample3,,,FALSE,
    ```

??? toggle "**Example use of Phandango coloring**"

    Data summary produced using the `phandango_coloring` option, visualized alongside Newick tree at <http://jameshadfield.github.io/phandango/#/main>

    !!! caption "Example phandango_coloring output"
        ![Phandango coloring example](../../assets/figures/example_phandango_coloring.png)

### Outputs

<div class="searchable-table" markdown="1">

| **Variable** | **Type** | **Description** |
|---|---|---|
| ksnp3_core_snp_matrix | File | The SNP matrix made with the core genome; formatted for Phandango if `phandango_coloring` input is `true` |
| ksnp3_core_snp_matrix_status | String | Will print either `The core SNP matrix was produced` OR `The core SNP matrix could not be produced` |
| ksnp3_core_snp_table | File | Formatted version of ksnp3_vcf_ref_genome file with only core SNPs, sorted by number of occurrences in the sample set |
| ksnp3_core_tree | File | The phylogenetic tree made with the core genome |
| ksnp3_docker | String | The docker image used |
| ksnp3_filtered_metadata | File | Optional output file with filtered metadata that is only produced if the optional `summarize_data` task is used. |
| ksnp3_ml_tree | File | Maximum likelihood tree that is only produced if `ksnp3_args` includes `"-ML"`  |
| ksnp3_nj_tree | File | Neighbor joining tree that is only produced if `ksnp3_args` includes `"-NJ"` |
| ksnp3_number_core_snps | String | Number of core SNPs in the sample set |
| ksnp3_number_snps | String | Number of SNPs in the sample set |
| ksnp3_pan_snp_matrix | File | The SNP matrix made with the pangenome; formatted for Phandango if `phandango_coloring` input is `true` |
| ksnp3_pan_tree | File | The phylogenetic tree made with the pangenome |
| ksnp3_snp_dists_version | String | The version of snp_dists used in the workflow |
| ksnp3_snps | File | File containing the set of SNPs used in the analysis. Required if more trees are to be appended to the existing one.  |
| ksnp3_summarized_data | File | CSV presence/absence matrix generated by the `summarize_data` task from the list of columns provided; formatted for Phandango if `phandango_coloring` input is `true` |
| ksnp3_vcf_ref_genome | File | A VCF file containing the variants detected in the core genome |
| ksnp3_vcf_ref_samplename | String | The name of the (user-supplied) sample used as the reference for calling SNPs. |
| ksnp3_vcf_snps_not_in_ref | File | A TSV file of the SNPs not present in the reference genome, but were identified by kSNP3. |
| ksnp3_wf_analysis_date | String | The date the workflow was run |
| ksnp3_wf_version | String | The version of the repository the workflow is hosted in |

</div>

## References

>Shea N Gardner, Tom Slezak, Barry G. Hall, kSNP3.0: SNP detection and phylogenetic analysis of genomes without genome alignment or reference genome, _Bioinformatics_, Volume 31, Issue 17, 1 September 2015, Pages 2877–2878, <https://doi.org/10.1093/bioinformatics/btv271>
<!-- -->
<https://github.com/tseemann/snp-dists>
