# Snippy_Tree

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Phylogenetic Construction](../../workflows_overview/workflows_type.md/#phylogenetic-construction) | [Bacteria](../../workflows_overview/workflows_kingdom.md/#bacteria) | PHB v2.1.0 | Yes; some optional features incompatible | Set-level |

## Snippy_Tree_PHB

`Snippy_Tree` is a workflow for generating high-quality bacterial phylogenies. It produces a phylogenetic tree and pairwise SNP-distance matrix, with the option to summarize additional metadata to visualize with the tree.

The tree produced by Snippy_Tree will always be a maximum-likelihood phylogeny using a reference-based alignment. There are key options for whether to:

- Generate a core-genome or whole-genome phylogeny (`core_genome`)
- Mask specified regions of the genome with a bed file (e.g. known repetitive regions for TB) (`bed_file`)
- Mask recombination (`use_gubbins`)
- Decide which nucleotide substitution model to use

### Inputs

`Snippy_Tree` is intended to be run after the `Snippy_Variants` workflow. It is a set-level workflow that takes in an array of directories generated by the `Snippy_Variants` workflow, which must be run for each sample that you wish to include in the phylogenetic tree. You should ensure that for all samples included in the phylogeny, `Snippy_Variants` has been run with identical inputs including the same reference genome. When running the `Snippy_Tree` workflow, you will need to provide the same reference genome that you used when running `Snippy_Variants`. `Snippy_Variants` and `Snippy_Tree` can both automatically be run by using the `Snippy_Streamline` workflow.

Sequencing data used in the Snippy_Tree workflow must:

- Be Illumina reads
- Be generated by unbiased whole genome shotgun sequencing
- Pass appropriate QC thresholds for the taxa to ensure that the reads represent reasonably complete genomes that are free of contamination from other taxa or cross-contamination of the same taxa.
- If masking recombination with `Gubbins`, input data should represent whole genomes from the same strain/lineage (e.g. MLST) that share a recent common ancestor.

!!! tip "Guidance for optional inputs"

    Several core and optional tasks can be used to generate the Snippy phylogenetic tree, making it highly flexible and suited to a wide range of datasets. You will need to decide which tasks to use depending on the genomes that you are analyzing. Some guidelines for the optional tasks to use for different genome types are provided below.
    
    ??? toggle "Default settings (suitable for most bacteria)"
    
        The default settings are as follows and are suitable for generating phylogenies for most bacteria
        
        - `core_genome` = true (creates core genome phylogeny)
        - `use_gubbins` = true (recombination masked)
        - nucleotide substitution model will be defined by IQTree's Model Finder
    
    ??? toggle "Phylogenies of _Mycobacterium tuberculosis_ complex"
    
        Phylogenies of MTBC are typically constructed
        
        - Using the H37Rv reference genome
            - `reference_genome_file` = gs://theiagen-public-files-rp/terra/theiaprok-files/Mtb_NC_000962.3.fasta
        - Masking repetitive regions of the genome (e.g. PE/PPE genes) that are often misaligned
            - `snippy_core_bed` = gs://theiagen-public-files/terra/theiaprok-files/Mtb_NC_000962.3.bed
        - Without masking recombination because TB can be considered non-recombinant
            - `use_gubbins` = false
        - Using the core genome
            - `core_genome` = true (as default)

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| snippy_tree_wf | **tree_name_updated** | String | Internal component, do not modify. Used for replacing spaces with underscores_ |  | Do not modify |
| snippy_tree_wf | **reference_genome_file** | File | Reference genome in FASTA or GENBANK format (must be the same reference used in Snippy_Variants workflow) |  | Required |
| snippy_tree_wf | **samplenames** | Array[String] | Samplenames for each input genome |  | Required |
| snippy_tree_wf | **snippy_variants_outdir_tarball** | Array[File] | Output from the Snippy_Variants workflow |  | Required |
| snippy_tree_wf | **tree_name** | String | String of your choice to prefix output files |  | Required |
| cg_reorder_matrix | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| cg_reorder_matrix | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| cg_reorder_matrix | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/mykrobe:0.12.1 | Optional |
| cg_reorder_matrix | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional |
| cg_snp_dists | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional |
| cg_snp_dists | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 50 | Optional |
| cg_snp_dists | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional |
| concatenate_variants | **docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1 | Optional |
| gubbins | **filter_percent** | Int | Maximum % gaps to include a sample in gubbins analysis and downstream analyses | 25 | Optional |
| gubbins | **iterations** | Int | Maximum number of trees to iteratively build to remove recombination | 5 | Optional |
| gubbins | **nuc_subst_model** | String | Nucleotide substitution model to use with Gubbins: "JC", "K2P", "HKY", "GTR", "GTRGAMMA" or "GTRCAT" (see <https://github.com/nickjcroucher/gubbins/blob/v3.3/docs/gubbins_manual.md#nucleotide-substitution-model-options>) | GTRGAMMA | Optional |
| gubbins | **tree_args** | String | Quoted string of further arguments passed to tree building algorithm |  | Optional |
| gubbins | **tree_builder** | String | Application to use for Gubbins tree building algorithm: "raxml", "raxmlng", "iqtree", "iqtree-fast", "fasttree", "hybrid" (fasttree is used for the first tree, and raxml is used for later iterations), "rapidnj" | raxml | Optional |
| iqtree2 | **alrt** | Int | Number of replicates to use for the SH-like approximate likelihood ratio test (Minimum recommended= 1000) | 1000 | Optional |
| shared_variants | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional |
| shared_variants | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| shared_variants | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16 | Optional |
| shared_variants | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| snippy_tree_wf | **call_shared_variants** | Boolean | When true, workflow generates table that combines variants across all samples and a table showing variants shared across samples | TRUE | Optional |
| snippy_tree_wf | **core_genome** | Boolean | When true, workflow generates core genome phylogeny; when false, whole genome is used | TRUE | Optional |
| snippy_tree_wf | **data_summary_column_names** | String | A comma-separated list of the column names from the sample-level data table for generating a data summary (presence/absence .csv matrix) |  | Optional |
| snippy_tree_wf | **data_summary_terra_project** | String | The billing project for your current workspace. This can be found after the "#workspaces/" section in the workspace's URL |  | Optional |
| snippy_tree_wf | **data_summary_terra_table** | String | The name of the sample-level Terra data table that will be used for generating a data summary |  | Optional |
| snippy_tree_wf | **data_summary_terra_workspace** | String | The name of the Terra workspace you are in. This can be found at the top of the webpage, or in the URL after the billing project. |  | Optional |
| snippy_tree_wf | **gubbins_cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional |
| snippy_tree_wf | **gubbins_docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/biocontainers/gubbins:3.3--py310pl5321h8472f5a_0 | Optional |
| snippy_tree_wf | **gubbins_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 32 | Optional |
| snippy_tree_wf | **gubbins_disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| snippy_tree_wf | **iqtree2_bootstraps** | String | Number of replicates for <http://www.iqtree.org/doc/Tutorial#assessing-branch-supports-with-ultrafast-bootstrap-approximation> (Minimum recommended= 1000) | 1000 | Optional |
| snippy_tree_wf | **iqtree2_cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional |
| snippy_tree_wf | **iqtree2_disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| snippy_tree_wf | **iqtree2_docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/iqtree2:2.1.2 | Optional |
| snippy_tree_wf | **iqtree2_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 32 | Optional |
| snippy_tree_wf | **iqtree2_model** | String | Nucelotide substitution model to use when generating the final tree with IQTree2. By default, IQtree runs its ModelFinder algorithm to identify the model it thinks best fits your dataset |  | Optional |
| snippy_tree_wf | **iqtree2_opts** | String | Additional options to pass to IQTree2 |  | Optional |
| snippy_tree_wf | **midpoint_root_tree** | Boolean | If true, midpoint root the final tree |  | Optional |
| snippy_tree_wf | **phandango_coloring** | Boolean | Boolean variable that tells the data summary task and the reorder matrix task to include a suffix that enables consistent coloring on Phandango; by default, this suffix is not added. To add this suffix set this variable to true. | FALSE | Optional |
| snippy_tree_wf | **snippy_core_bed** | File | Bed file with locations to be masked from the core genome alignment |  | Optional |
| snippy_tree_wf | **snippy_core_cpu** | Int | Number of CPUs to allocate to the task | 8 | Optional |
| snippy_tree_wf | **snippy_core_disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| snippy_tree_wf | **snippy_core_docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/snippy:4.6.0 | Optional |
| snippy_tree_wf | **snippy_core_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 16 | Optional |
| snippy_tree_wf | **snp_dists_docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/snp-dists:0.8.2 | Optional |
| snippy_tree_wf | **snp_sites_cpus** | Int | CPUs to allocate to SNP-sites | 1 | Optional |
| snippy_tree_wf | **snp_sites_disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| snippy_tree_wf | **snp_sites_docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/snp-sites:2.5.1 | Optional |
| snippy_tree_wf | **snp_sites_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 4 | Optional |
| snippy_tree_wf | **use_gubbins** | Boolean | When "true", workflow removed recombination with gubbins tasks; when "false", gubbins is not used | true | Optional |
| summarize_data | **cpu** | Int | Number of CPUs to allocate to the task | 8 | Optional |
| summarize_data | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| summarize_data | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16 | Optional |
| summarize_data | **id_column_name** | String | Name of the column in the input table that contains the sample IDs, if different from default |  | Optional |
| summarize_data | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 1 | Optional |
| version_capture | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional |
| wg_reorder_matrix | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| wg_reorder_matrix | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| wg_reorder_matrix | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/mykrobe:0.12.1 | Optional |
| wg_reorder_matrix | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional |
| wg_snp_dists | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional |
| wg_snp_dists | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 50 | Optional |
| wg_snp_dists | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional |

### Workflow Tasks

??? task "Snippy"

    ##### Snippy {#snippy_task}

    Snippy is a pipeline for calling SNPs and INDELs in haploid genomes. Before running `Snippy_Tree`, you must run `Snippy_Variants`, another workflow that uses the Snippy tool to align reads against a reference genome for individual samples. In `Snippy_Tree`, the snippy tool is used again to generate a whole-genome multiple sequence alignment (fasta file) of reads from all the samples we'd like in our tree. 

    When generating the multiple sequence alignment, a bed file can be provided by users to mask certain areas of the genome in the alignment. This is particularly relevant for masking known repetitive regions in _Mycobacterium tuberculosis_  genomes, or masking known regions containing phage sequences.

    !!! info "Why do I see `snippy_core` in Terra?"
        In Terra, this task is named "snippy_core" after the name of the command in the original Snippy tool. Despite the name, this command is NOT being used to make a core genome, but instead a multiple sequence alignment of the whole genome (without any sections masked using a bed file).
        
    !!! techdetails "Snippy Technical Details"
    
        |  | Links |
        | --- | --- |
        | Task | [task_snippy_core.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/utilities/task_snippy_core.wdl) |
        | Default software version | v4.6.0 (us-docker.pkg.dev/general-theiagen/staphb/snippy:4.6.0) |
        | Software Source Code | [Snippy on GitHub](https://github.com/tseemann/snippy) |
        | Software Documentation | [Snippy on GitHub](https://github.com/tseemann/snippy) |

??? task "Gubbins (optional)"

    ##### Gubbins (optional) {#gubbins_task}

    !!! info "Optional"
        Gubbins is used when `use_gubbins` is set to `true` (default=true).

    **G**enealogies **U**nbiased **B**y recom**B**inations **I**n **N**ucleotide **S**equences (Gubbins) identifies and masks genomic regions that are predicted to have arisen via recombination. It works by iteratively identifying loci containing elevated densities of SNPs and constructing phylogenies based on the putative single nucleotide variants outside these regions (for more details, see [here](https://github.com/nickjcroucher/gubbins/blob/v3.3/docs/gubbins_manual.md#description-of-the-algorithm)). By default, these phylogenies are constructed using RaxML and a GTR-GAMMA nucleotide substitution model, which will be the most suitable model for most bacterial phylogenetics, though this can be modified with the `tree_builder` and `nuc_subst_model` inputs.

    Gubbins is the industry standard for masking recombination from bacterial genomes when building phylogenies, but limitations to recombination removal exist. Gubbins cannot distinguish recombination from high densities of SNPs that may result from assembly or alignment errors, mutational hotspots, or regions of the genome with relaxed selection. The tool is also intended only to find recombinant regions that are short relative to the length of the genome, so large regions of recombination may not be masked. These factors should be considered when interpreting resulting phylogenetic trees, but overwhelmingly Gubbins improves our ability to understand ancestral relationships between bacterial genomes.

    There are few optional inputs for Gubbins that can be modified by the user:

    - `iterations`: Gubbins works by iteratively identifying loci containing elevated densities of SNPs, while constructing phylogenies based on the putative single nucleotide variants outside these regions. It may take many iterations for Gubbins to converge on an alignment that it considers free of recombination, especially for phylogenies that contain large numbers of genomes. By default, Gubbins is limited to 5 iterations though this may be increased by the user with the `iterations`optional input (incurring increased computing time and cost, and possibly requiring increased memory allocation).
    - `nuc_subst_model`, `tree_builder` and `tree_args`:  When Gubbins constructs phylogenies, it can use a number of phylogenetic inference tools, each with [different nucleotide substitution models](https://github.com/nickjcroucher/gubbins/blob/master/docs/gubbins_manual.md#nucleotide-substitution-model-options) and [tree-building models](https://github.com/nickjcroucher/gubbins/blob/master/docs/gubbins_manual.md#tree-building-options). By default, the `Snippy_Tree` workflow uses a GTRGAMMA substitution model and RaxML for tree building (typically suitable for bacterial genomes), but these can be modified by the user depending on the genome sequences being used with the `nuc_subst_model` and `tree_builder` optional inputs, respectively. The nucleotide substitution models that are available depend on the tree building algorithm being used (see [here](https://github.com/nickjcroucher/gubbins/blob/v3.3/docs/gubbins_manual.md#nucleotide-substitution-model-options)). Additional options for generating the phylogenetic trees in Gubbins can be specified with the `tree_args` optional input, providing an input string that is consistent with the option formats of the Gubbins command.
    - `filter_percent`: By default, Gubbins removes genomes from the multiple sequence alignment if  more than 25 % of the genome is represented by gaps. The percentage of gaps can be modified by the user using the `filter_percent` optional input.

    !!! techdetails "Gubbins Technical Details"
                
        |  | Links |
        | --- | --- |
        | Task | [task_gubbins.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/task_gubbins.wdl) |
        | Software Source Code | [Gubbins on GitHub](https://github.com/nickjcroucher/gubbins) |
        | Software Documentation | [Gubbins v3.3 manual](https://github.com/nickjcroucher/gubbins/blob/v3.3/docs/gubbins_manual.md) |
        | Original Publication(s) | [Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences using Gubbins](https://academic.oup.com/nar/article/43/3/e15/2410982) |
        | Default software version | us-docker.pkg.dev/general-theiagen/biocontainers/gubbins:3.3--py310pl5321h8472f5a_0 |

??? task "SNP-sites (optional)"

    ##### SNP-sites (optional) {#snp_sites_task}

    !!! tip "Turn on SNP-Sites with `core_genome`"
        SNP-sites runs when the `core_genome` option is set to true.

    SNP-sites is used to filter out invariant sites in the whole-genome alignment, thereby creating a core genome alignment for phylogenetic inference. The output is a fasta file containing the core genome of each sample only. If Gubbins has been used, this output fasta will not contain any sites that are predicted to have arisen via recombination.

    !!! techdetails "SNP-sites technical details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_snp_sites.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/utilities/task_snp_sites.wdl) |
        | Default software version | 2.5.1 (us-docker.pkg.dev/general-theiagen/biocontainers/snp-sites:2.5.1--hed695b0_0) |
        | Software Source Code | [SNP-sites on GitHub](https://github.com/sanger-pathogens/snp-sites) |
        | Software Documentation | [SNP-sites on GitHub](https://github.com/sanger-pathogens/snp-sites) |
        | Original Publication(s) | [SNP-sites: rapid efficient extraction of SNPs from multi-FASTA alignments](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000056) |

??? task "IQTree2"

    ##### IQTree2 {#iqtree2_task}

    IQTree2 is used to build the final phylogeny. It uses the alignment generated in the previous steps of the workflow. The contents of this alignment will depend on whether any sites were masked with recombination.

    The phylogeny is generated using the maximum-likelihood method and a specified nucleotide substitution model. By default, the Snippy_Tree workflow will run Model Finder to determine the most appropriate nucleotide substitution model for your data, but you may specify the nucleotide substitution model yourself using the `iqtree2_model` optional input (see [here](http://www.iqtree.org/doc/Substitution-Models) for available models).

    IQTree will perform assessments of the tree using the Shimodaira–Hasegawa approximate likelihood-ratio test ([SH-aLRT test](https://academic.oup.com/sysbio/article/59/3/307/1702850?login=false)), and ultrafast bootstrapping with [UFBoot2](https://academic.oup.com/mbe/article/35/2/518/4565479), a quicker but less biased alternative to standard bootstrapping. A clade should not typically be trusted if it has less than 80% support from the SH-aLRT test and less than 95% support with ultrafast bootstrapping.

    !!! tip "Nucleotide substitution model"
        When `core_genome`= `true`, the default nucleotide substitution model is set to the General Time Reverside model with Gamma distribution (GTR+G). 
        
        When the user sets `core_genome`= `false`, the default nucleotide substitution model is set to the General Time Reversible model with invariant sites and Gamma distribution (`GTR+I+G`).
                
    !!! techdetails "IQTree2 technical details"
                    
        |  | Links |
        | --- | --- |
        | Task | [task_iqtree2.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/task_iqtree2.wdl) |
        | Software Source Code | [IQ-TREE on GitHub](https://github.com/iqtree/iqtree2) |
        | Software Documentation | [IQTree documentation](http://www.iqtree.org/doc/) for the latest version (not necessarily the version used in this workflow) |
        | Original Publication(s) | [IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era](https://academic.oup.com/mbe/article/37/5/1530/5721363) |
        | Publication for the SH-alRT test | [New Algorithms and Methods to Estimate Maximum-Likelihood Phylogenies: Assessing the Performance of PhyML 3.0](https://academic.oup.com/sysbio/article/59/3/307/1702850?login=false) |
        | Publication for ultrafast bootstrapping integration to IQTree | [Ultrafast Approximation for Phylogenetic Bootstrap](https://academic.oup.com/mbe/article/30/5/1188/997508?login=false); [UFBoot2: Improving the Ultrafast Bootstrap Approximation](https://academic.oup.com/mbe/article/35/2/518/4565479?login=false) |
        | Publication for ModelFinder  | [ModelFinder: fast model selection for accurate phylogenetic estimates](https://www.nature.com/articles/nmeth.4285) |

??? task "SNP-dists"

    ##### SNP-dists {#snp_dists_task}

    `SNP-dists` computes pairwise SNP distances between genomes. It takes the same alignment of genomes used to generate your phylogenetic tree and produces a matrix of pairwise SNP distances between sequences. This means that if you generated pairwise core-genome phylogeny, the output will consist of pairwise core-genome SNP (cgSNP) distances. Otherwise, these will be whole-genome SNP distances. Regardless of whether core-genome or whole-genome SNPs, this SNP distance matrix will exclude all SNPs in masked regions (i.e. masked with a bed file or gubbins). 

    The SNP-distance output can be visualized using software such as [Phandango](http://jameshadfield.github.io/phandango/#/main) to explore the relationships between the genomic sequences. The task adds a Phandango coloring tag (:c1) to the column names in the output matrix to ensure that all columns are colored with the same color scheme throughout.

    !!! techdetails "SNP-dists Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_snp_dists.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/utilities/task_snp_dists.wdl) |
        | Default software version | 0.8.2 (us-docker.pkg.dev/general-theiagen/staphb/snp-dists:0.8.2) |
        | Software Source Code | [SNP-dists on GitHub](https://github.com/tseemann/snp-dists) |
        | Software Documentation | [SNP-dists on GitHub](https://github.com/tseemann/snp-dists) |
        | Original Publication(s) | Not known to be published |

??? task "Data summary (optional)"

    ##### Data Summary (optional) {#data_summary_task}

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

    !!! techdetails "Data summary technical details"

        |  | Links |
        | --- | --- |
        | Task | [task_summarize_data.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/task_summarize_data.wdl) |

??? task "Concatenate Variants (optional)"

    ##### Concatenate Variants (optional) {#concatenate_variants_task}

    The `cat_variants` task concatenates variant data from multiple samples into a single file `concatenated_variants`. It is very similar to the `cat_files` task, but also adds a column to the output file that indicates the sample associated with each row of data.

    The `concatenated_variants` file will be in the following format:

    | samplename | CHROM | POS | TYPE | REF | ALT | EVIDENCE | FTYPE | STRAND | NT_POS | AA_POS | EFFECT | LOCUS_TAG | GENE | PRODUCT |
    | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
    | sample1 | PEKT02000007 | 5224 | snp | C | G | G:21 C:0 |  |  |  |  |  |  |  |  |
    | sample2 | PEKT02000007 | 34112 | snp | C | G | G:32 C:0 | CDS | + | 153/1620 | 51/539 | missense_variant c.153C>G p.His51Gln | B9J08_002604 | hypothetical protein |  |
    | sample3 | PEKT02000007 | 34487 | snp | T | A | A:41 T:0 | CDS | + | 528/1620 | 176/539 | missense_variant c.528T>A p.Asn176Lys | B9J08_002604 | hypothetical protein |  |

    !!! techdetails "Technical Details"
    
        |  | Links |
        | --- | --- |
        | Task | /tasks/utilities/file_handling/task_cat_files.wdl |
        | Software Source Code | [task_cat_files.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/file_handling/task_cat_files.wdl) |

??? task "Shared Variants Task (Optional)"

    ##### Shared Variants (optional) {#shared_variants_task}

    The `shared_variants` task takes in the `concatenated_variants` output from the `cat_variants` task and reshapes the data so that variants are rows and samples are columns. For each variant, samples where the variant was detected are populated with a "1" and samples were **either the variant was not detected or there was insufficient coverage to call variants** are populated with a "0". The resulting table is available as the `shared_variants_table` output.

    The `shared_variants_table` file will be in the following format:

    | CHROM | POS | TYPE | REF | ALT | FTYPE | STRAND | NT_POS | AA_POS | EFFECT | LOCUS_TAG | GENE | PRODUCT | sample1 | sample2 | sample3 |
    | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
    | PEKT02000007 | 2693938 | snp | T | C | CDS | - | 1008/3000 | 336/999 | synonymous_variant c.1008A>G p.Lys336Lys | B9J08_003879 | NA | chitin synthase 1 | 1 | 1 | 0 |
    | PEKT02000007 | 2529234 | snp | G | C | CDS | + | 282/336 | 94/111 | missense_variant c.282G>C p.Lys94Asn | B9J08_003804 | NA | cytochrome c | 1 | 1 | 1 |
    | PEKT02000002 | 1043926 | snp | A | G | CDS | - | 542/1464 | 181/487 | missense_variant c.542T>C p.Ile181Thr | B9J08_000976 | NA | dihydrolipoyl dehydrogenase | 1 | 1 | 0 |
    
    !!! techdetails "Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | task_shared_variants.wdl |
        | Software Source Code | [task_shared_variants.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/utilities/task_shared_variants.wdl) |

### Outputs

| **Variable** | **Type** | **Description** |
|---|---|---|
| snippy_cg_snp_matrix | File | CSV file of core genome pairwise SNP distances between samples, calculated from the final alignment  |
| snippy_concatenated_variants | File | Concatenated snippy_results file across all samples in the set |
| snippy_filtered_metadata | File | TSV recording the columns of the Terra data table that were used in the summarize_data task |
| snippy_final_alignment | File | Final alignment (FASTA file) used to generate the tree (either after snippy alignment, gubbins recombination removal, and/or core site selection with SNP-sites) |
| snippy_final_tree | File | Newick tree produced from the final alignment. Depending on user input for core_genome, the tree could be a core genome tree (default when core_genome is true) or whole genome tree (if core_genome is false) |
| snippy_gubbins_branch_stats | File | CSV file showing https://github.com/nickjcroucher/gubbins/blob/master/docs/gubbins_manual.md#output-statistics for each branch of the tree   |
| snippy_gubbins_docker | String | Docker file used for running Gubbins |
| snippy_gubbins_recombination_gff | File | Recombination statistics in GFF format; these can be viewed in Phandango against the phylogenetic tree |
| snippy_gubbins_version | String | Gubbins version used |
| snippy_iqtree2_docker | String | Docker file used for running IQTree2 |
| snippy_iqtree2_model_used | String | Nucleotide substitution model used by IQTree2 |
| snippy_iqtree2_version | String |  IQTree2 version used |
| snippy_msa_snps_summary | File | TXT file containing summary statistics for each alignment of each input genome against the reference. This indicates how good the alignment is. Pay particular attention to # unaligned sites, and heterogeneous positions. |
| snippy_ref | File | Reference genome (FASTA or GenBank file) used for generating phylogeny |
| snippy_shared_snp_table | File | Table illustrating variants shared among samples |
| snippy_snp_dists_docker | String | Docker file used for running SNP-dists |
| snippy_snp_dists_version | String | SNP-dists version used |
| snippy_snp_sites_docker | String | Docker file used for running SNP-sites  |
| snippy_snp_sites_version | String | SNP-sites version used |
| snippy_summarized_data | File | CSV presence/absence matrix generated by the summarize_data task from the list of columns provided; formatted for Phandango if phandango_coloring input is true |
| snippy_tree_analysis_date | String | Date of workflow run |
| snippy_tree_snippy_docker | String | Docker file used for running Snippy |
| snippy_tree_snippy_version | String | Snippy version used |
| snippy_tree_version | String | Version of Snippy_Tree workflow |
| snippy_wg_snp_matrix | File | CSV file of whole genome pairwise SNP distances between samples, calculated from the final alignment |

## References

> **Gubbins:** Croucher, Nicholas J., Andrew J. Page, Thomas R. Connor, Aidan J. Delaney, Jacqueline A. Keane, Stephen D. Bentley, Julian Parkhill, and Simon R. Harris. 2015. "Rapid Phylogenetic Analysis of Large Samples of Recombinant Bacterial Whole Genome Sequences Using Gubbins." Nucleic Acids Research 43 (3): e15.
<!-- -->
> **SNP-sites:** Page, Andrew J., Ben Taylor, Aidan J. Delaney, Jorge Soares, Torsten Seemann, Jacqueline A. Keane, and Simon R. Harris. 2016. "SNP-Sites: Rapid Efficient Extraction of SNPs from Multi-FASTA Alignments." Microbial Genomics 2 (4): e000056.
<!-- -->
> **IQTree:** Nguyen, Lam-Tung, Heiko A. Schmidt, Arndt von Haeseler, and Bui Quang Minh. 2015. "IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies." Molecular Biology and Evolution 32 (1): 268–74.