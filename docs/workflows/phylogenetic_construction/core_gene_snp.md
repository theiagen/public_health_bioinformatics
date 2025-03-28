# Core_Gene_SNP

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Phylogenetic Construction](../../workflows_overview/workflows_type.md/#phylogenetic-construction) | [Bacteria](../../workflows_overview/workflows_kingdom.md/#bacteria) | PHB v3.0.0 | Yes, some optional features incompatible | Set-level |

## Core_Gene_SNP_PHB

!!! caption "Core Gene SNP Workflow Diagram"
    ![Core Gene SNP Workflow Diagram](../../assets/figures/Core_Gene_SNP.png){width:45%}

The Core_Gene_SNP workflow is intended for pangenome analysis, core gene alignment, and phylogenetic analysis. The workflow takes in gene sequence data in GFF3 format from a set of samples. It first produces a pangenome summary using [`Pirate`](https://github.com/SionBayliss/PIRATE), which clusters genes within the sample set into orthologous gene families. By default, the workflow also instructs `Pirate` to produce both core gene and pangenome alignments. The workflow subsequently triggers the generation of a phylogenetic tree and SNP distance matrix from the core gene alignment using [`iqtree`](https://github.com/iqtree/iqtree2/tree/v1.6.7) and [`snp-dists`](https://github.com/tseemann/snp-dists), respectively. Optionally, the workflow will also run this analysis using the pangenome alignment. This workflow also features an optional module, `summarize_data`, that creates a presence/absence matrix for the analyzed samples from a list of indicated columns (such as AMR genes, etc.) that can be used in Phandango.

!!! info "Default Parameters"
    Please note that while default parameters for pangenome construction and phylogenetic tree generation are provided, **these default parameters may not suit every dataset and have not been validated against known phylogenies**. Users should take care to select the parameters that are most appropriate for their dataset. Please reach out to [support@theiagen.com](mailto:support@theiagen.com) or one of the other resources listed at the bottom of this page if you would like assistance with this task.

### Inputs

For further detail regarding Pirate options, please see [PIRATE's documentation](https://github.com/SionBayliss/PIRATE). For further detail regarding IQ-TREE options, please see `http://www.iqtree.org/doc/Command-Reference`.

This workflow runs on the set level.

<div class="searchable-table" markdown="1">

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| core_gene_snp_workflow | **cluster_name** | String | Name of sample set | | Required |
| core_gene_snp_workflow | **gff3** | Array[File] | Array of gff3 files to include in analysis, output gff files from both prokka and bakta using TheiaProk workflows are compatible | | Required |
| core_gene_snp_workflow | **midpoint_root_tree** | Boolean | Boolean variable that will instruct the workflow to reroot the tree at the midpoint | FALSE | Optional |
| core_gene_snp_workflow | **phandango_coloring** | Boolean | Boolean variable that tells the data summary task and the reorder matrix task to include a suffix that enables consistent coloring on Phandango; by default, this suffix is not added. To add this suffix set this variable to true. | FALSE | Optional |
| core_gene_snp_workflow | **data_summary_terra_table** | String | The name of the Terra data table that you want data pulled from |  | Optional |
| core_gene_snp_workflow | **data_summary_column_names** | String | A comma-delimited list of columns in the origin data table that contains contain that you would like a presence/absence .csv matrix generated for |  | Optional |
| core_gene_snp_workflow | **core_tree** | Boolean | Boolean variable that instructs the workflow to create a phylogenetic tree and SNP distance matrix from the core gene alignment. Align must also be set to true. | TRUE | Optional |
| core_gene_snp_workflow | **pan_tree** | Boolean | Boolean variable that instructs the workflow to create a phylogenetic tree and SNP distance matrix from the pangenome alignment. Align must also be set to true. | FALSE | Optional |
| core_gene_snp_workflow | **data_summary_terra_workspace** | String | The name of the current Terra workspace you are in; this can be found at the top of the webpage, or in the URL after the billing project. |  | Optional |
| core_gene_snp_workflow | **align** | Boolean | Boolean variable that instructs the workflow to generate core and pangenome alignments if "true". If "false", the workflow will produce only a pangenome summary. | TRUE | Optional |
| core_gene_snp_workflow | **data_summary_terra_project** | String | The billing project for the current workspace; can be found after the "#workspaces/" section in the workflow's URL |  | Optional |
| core_gene_snp_workflow | **sample_names** | Array[String] | Array of sample_ids from the data table used |  | Optional |
| core_iqtree | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 32 | Optional |
| core_iqtree | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| core_iqtree | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional |
| core_iqtree | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/iqtree:1.6.7 | Optional |
| core_iqtree | **iqtree_model** | String | Substitution model, frequency type (optional) and rate heterogeneity type (optional) used by IQ-TREE. This string follows the IQ-TREE "-m" option. For comparison to other tools use HKY for Bactopia, GTR+F+I for Grandeur, GTR+G4 for Nullarbor, GTR+G for Dryad | GTR+I+G | Optional |
| core_iqtree | **iqtree_opts** | String | Additional options for IQ-TREE, see <http://www.iqtree.org/doc/Command-Reference> |  | Optional |
| core_iqtree | **iqtree_bootstraps** | String | Number of ultrafast bootstrap replicates. Follows IQ-TREE "-bb" option. | 1000 | Optional |
| core_iqtree | **alrt** | String | Number of replicates to perform SH-like approximate likelihood ratio test (SH-aLRT). Follows IQ-TREE "-alrt" option | 1000 | Optional |
| core_reorder_matrix | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| core_reorder_matrix | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| core_reorder_matrix | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/mykrobe:0.12.1 | Optional |
| core_reorder_matrix | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional |
| core_snp_dists | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional |
| core_snp_dists | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/snp-dists:0.8.2 | Optional |
| core_snp_dists | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 50 | Optional |
| core_snp_dists | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional |
| pan_iqtree | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 32 | Optional |
| pan_iqtree | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| pan_iqtree | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional |
| pan_iqtree | **alrt** | String | Number of replicates to perform SH-like approximate likelihood ratio test (SH-aLRT). Follows IQ-TREE "-alrt" option | 1000 | Optional |
| pan_iqtree | **iqtree_model** | String | Substitution model, frequency type (optional) and rate heterogeneity type (optional) used by IQ-TREE. This string follows the IQ-TREE "-m" option. For comparison to other tools use HKY for Bactopia, GTR+F+I for Grandeur, GTR+G4 for Nullarbor, GTR+G for Dryad | GTR+I+G | Optional |
| pan_iqtree | **iqtree_bootstraps** | String | Number of ultrafast bootstrap replicates. Follows IQ-TREE "-bb" option. | 1000 | Optional |
| pan_iqtree | **iqtree_opts** | String | Additional options for IQ-TREE, see <http://www.iqtree.org/doc/Command-Reference> |  | Optional |
| pan_iqtree | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/iqtree:1.6.7 | Optional |
| pan_reorder_matrix | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| pan_reorder_matrix | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| pan_reorder_matrix | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/mykrobe:0.12.1 | Optional |
| pan_reorder_matrix | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional |
| pan_snp_dists | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional |
| pan_snp_dists | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 50 | Optional |
| pan_snp_dists | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/snp-dists:0.8.2 | Optional |
| pan_snp_dists | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional |
| pirate | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| pirate | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional |
| pirate | **nucl** | Boolean | Boolean variable that instructs pirate to create a pangenome on CDS features using nucleotide identity, rather than amino acid identity, if true.  | FALSE | Optional |
| pirate | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 32 | Optional |
| pirate | **panopt** | String | Additional arguments for Pirate |  | Optional |
| pirate | **docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/biocontainers/pirate:1.0.5--hdfd78af_0 | Optional |
| pirate | **features** | String | Features to use for pangenome construction [default: CDS] | CDS | Optional |
| pirate | **steps** | String | Identity thresholds to use for pangenome construction | 50,60,70,80,90,95,98 | Optional |
| summarize_data | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| summarize_data | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16 | Optional |
| summarize_data | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 1 | Optional |
| summarize_data | **id_column_name** | String | Use in the case your sample IDs are not in the table ID column | 1 | Optional |
| summarize_data | **cpu** | Int | Number of CPUs to allocate to the task | 8 | Optional |
| version_capture | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional |

</div>

### Workflow Tasks

By default, the Core_Gene_SNP workflow will begin by analyzing the input sample set using [PIRATE](https://github.com/SionBayliss/PIRATE). Pirate takes in GFF3 files and classifies the genes into gene families by sequence identity, outputting a pangenome summary file. The workflow will instruct Pirate to create core gene and pangenome alignments using this gene family data. Setting the "align" input variable to false will turn off this behavior, and the workflow will output only the pangenome summary. The workflow will then use the core gene alignment from `Pirate` to infer a phylogenetic tree using `IQ-TREE`. It will also produce an SNP distance matrix from this alignment using [snp-dists](https://github.com/tseemann/snp-dists). This behavior can be turned off by setting the `core_tree` input variable to false. The workflow will not create a pangenome tree or SNP-matrix by default, but this behavior can be turned on by setting the `pan_tree` input variable to true.

The optional `summarize_data` task performs the following only if all of the `data_summary_*` and `sample_names` optional variables are filled out:

1. Digests a _comma-separated_  list of column names, such as `"amrfinderplus_virulence_genes,amrfinderplus_stress_genes"`, etc. that can be found within the origin Terra data table.
2. It will then parse through those column contents and extract each value; for example, if the `amrfinder_amr_genes` column for a sample contains these values: `"aph(3')-IIIa,tet(O),blaOXA-193"`, the `summarize_data` task will check each sample in the set to see if they also have those AMR genes detected.
3. Outputs a .csv file that indicates presence (TRUE) or absence (empty) for each item in those columns; that is, it will check each sample in the set against the detected items in each column to see if that value was also detected.

By default, this task appends a Phandango coloring tag to color all items from the same column the same; this can be turned off by setting the optional `phandango_coloring` variable to `false`.

### Outputs

<div class="searchable-table" markdown="1">

| **Variable** | **Type** | **Description** |
|---|---|---|
| core_gene_snp_wf_analysis_date | String | Date of analysis using Core_Gene_SNP workflow |
| core_gene_snp_wf_version | String | Version of PHBG used for analysis |
| pirate_core_alignment_fasta | File | Nucleotide alignments of the core genes as created using MAFFT within Pirate. Loci are ordered according to the gene_families.ordered file. |
| pirate_core_alignment_gff | File | Annotation data for the gene family within the corresponding fasta file |
| pirate_core_snp_matrix | File | SNP distance matrix created from the core gene alignment |
| pirate_docker_image | String | Pirate docker image used |
| pirate_gene_families_ordered | File | Summary of all gene families, as estimated by Pirate |
| pirate_iqtree_core_tree | File | Phylogenetic tree produced by IQ-TREE from the core gene alignment |
| pirate_iqtree_pan_tree | File | Phylogenetic tree produced by IQ-TREE from the pangenome alignment |
| pirate_iqtree_version | String | IQ-TREE version used |
| pirate_pan_alignment_fasta | File | Nucleotide alignments of the pangenome by gene as created using MAFFT within Pirate. Loci are ordered according to the gene_families.ordered file. |
| pirate_pan_alignment_gff | File | Annotation data for the gene family within the corresponding fasta file |
| pirate_pan_snp_matrix | File | SNP distance matrix created from the pangenome alignment |
| pirate_pangenome_summary | File | Summary of the number and frequency of genes in the pangenome, as estimated by Pirate |
| pirate_presence_absence_csv | File | A file generated by Pirate that allows many post-alignment tools created for Roary to be used on the output from Pirate |
| pirate_snp_dists_version | String | Version of snp-dists used  |
| pirate_summarized_data | File | The presence/absence matrix generated by the summarize_data task from the list of columns provided |

</div>

## References

>Sion C Bayliss, Harry A Thorpe, Nicola M Coyle, Samuel K Sheppard, Edward J Feil, PIRATE: A fast and scalable pangenomics toolbox for clustering diverged orthologues in bacteria, _GigaScience_, Volume 8, Issue 10, October 2019, giz119, <https://doi.org/10.1093/gigascience/giz119>
<!-- -->
> Lam-Tung Nguyen, Heiko A. Schmidt, Arndt von Haeseler, Bui Quang Minh, IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies, _Molecular Biology and Evolution_, Volume 32, Issue 1, January 2015, Pages 268–274, <https://doi.org/10.1093/molbev/msu300>
<!-- -->
> <https://github.com/tseemann/snp-dists>
