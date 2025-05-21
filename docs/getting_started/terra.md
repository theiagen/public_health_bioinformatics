# Getting Started with Terra

!!! dna "Our Approach"
    Theiagen’s approach to genomic analysis in public health typically uses the [Terra](https://terra.bio/) platform to run workflows that undertake bioinformatic analysis, then uses other platforms for visualization of the resulting data. This is described in more depth in our paper [_Accelerating bioinformatics implementation in public health_](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001051), and the application of this approach for genomic surveillance of SARS-CoV-2 in California is described in the paper [_Pathogen genomics in public health laboratories: successes, challenges, and lessons learned from California’s SARS-CoV-2 Whole-Genome Sequencing Initiative, California COVIDNet_](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001027).

!!! tip ""
    **When undertaking genomic analysis using Terra and other data visualization platforms, it is essential to consider the necessary and appropriate workflows and resources for your analysis. To help you make these choices, take a look at the relationship between the most commonly used Theiagen workflows, and the descriptions of the major stages in genomic data analysis below.**

    !!! caption "Analysis Approaches for Genomic Data"
        ![The relationship between the various PHB workflows](../assets/figures/Workflow_Relationships.png#only-light){data-description="This diagram shows the Theiagen workflows (green boxes) available for analysis of genomic data in public health and the workflows that may be used consecutively (arrows). The blue boxes describe the major functions that these workflows undertake. The yellow boxes show functions that may be undertaken independently of workflows on Terra."}
        ![The relationship between the various PHB workflows](../assets/figures/Workflow_Relationships_dark.png#only-dark){data-description="This diagram shows the Theiagen workflows (green boxes) available for analysis of genomic data in public health and the workflows that may be used consecutively (arrows). The blue boxes describe the major functions that these workflows undertake. The yellow boxes show functions that may be undertaken independently of workflows on Terra."}

        This diagram shows the Theiagen workflows (green boxes) available for analysis of genomic data in public health and the workflows that may be used consecutively (arrows). The blue boxes describe the major functions that these workflows undertake. The yellow boxes show functions that may be undertaken independently of workflows on Terra.

## Data Import to Terra

To start using Terra for data analysis, you will first need to import your data into your workspace. There are multiple ways to do this:

- **Using Terra’s native features to upload data from your local computer or link to data that’s already in a Google bucket**
- Data import workflows
    - Using the [SRA_Fetch](../workflows/data_import/sra_fetch.md) workflow to import publicly available data from any repository in the [INSDC](https://www.insdc.org/) (including with [SRA](https://www.ncbi.nlm.nih.gov/sra), [ENA](https://www.ebi.ac.uk/ena/browser/home) and [DRA](https://www.ddbj.nig.ac.jp/dra/index-e.html))
    - Using the [Assembly_Fetch ](../workflows/data_import/assembly_fetch.md) workflow to import publicly available genome assemblies from [NCBI](https://www.ncbi.nlm.nih.gov/datasets/)
    - Using the [BaseSpace_Fetch](../workflows/data_import/basespace_fetch.md) workflow to import data from your [Illumina BaseSpace](https://basespace.illumina.com/) account
    - Using the [Create_Terra_Table](../workflows/data_import/create_terra_table.md) workflow to help create your data table after manual upload to your Terra workspace (or a Google Cloud Storage Bucket)

!!! example "SOPs for ==importing data into a Terra workspace=="

    | SOP | SOP Version | PHB Version Compatibility | 
    |---|---|---|
    | [Uploading Data, Creating Metadata Tables and TSV files, and Importing Workflows](../assets/sops/TG-TER-03_GettingStartedInTerra_v3.pdf) | v3 | v1.3.0, v2+ |
    | [Linking BaseSpace and Importing BaseSpace Reads to Terra](../assets/sops/TG-TER-04_BaseSpaceFetch_v2.pdf) | v3 | v1.3.0, v2+ |

## Genome assembly, QC, and characterization

### TheiaX workflows

The TheiaX workflows are used for genome assembly, quality control, and characterization. The [TheiaCoV Workflow Series](../workflows/genomic_characterization/theiacov.md), [TheiaProk Workflow Series](../workflows/genomic_characterization/theiaprok.md), and [TheiaEuk Workflow Series](../workflows/genomic_characterization/theiaeuk.md) workflows are intended for viral, bacterial, and fungal pathogens, respectively. [TheiaMeta Workflow Series](../workflows/genomic_characterization/theiameta.md)  is intended for the analysis of a single taxon from metagenomic data.

!!! example "SOPs for the ==TheiaX workflows=="

    ??? toggle "For analyzing ==SARS-CoV-2=="
        | SOP | SOP Version | PHB Version Compatibility | 
        |---|---|---|
        | [Analyze SARS-COV-2 using TheiaCoV_Illumina_PE_PHB](../assets/sops/TG-SC2-PE_SC2_TheiaCoV_IlluminaPE_v3.pdf) | v3 | v2+ |
        | [Analyze SARS-COV-2 using TheiaCoV_Illumina_SE_PHB](../assets/sops/TG-SC2-SE_SC2_TheiaCoV_IlluminaSE_v3.pdf) | v3 | v2+ |
        | [Analyze SARS-COV-2 using TheiaCoV_ClearLabs](../assets/sops/TG-SC2-CL_SC2_ClearLabs_v3.pdf) | v3 | v2+ |
        | [Analyze SARS-COV-2 using TheiaCoV_ONT](../assets/sops/TG-SC2-ONT_Analyzing_SC2_Using_TheiaCov_ONT_PHB_v2.pdf) | v2 | v1.x+ |
        | [Analyzing SARS-CoV-2 using TheiaCoV_FASTA](../assets/sops/TG-SC2-FST_Analyzing_SC2_Using_TheiaCoV_FASTA_PHB_v2.pdf) | v2 | v1.x+ |

    ??? toggle "For analyzing ==influenza=="
        | SOP | SOP Version | PHB Version Compatibility | 
        |---|---|---|
        | [Analyzing Flu Data in Terra using TheiaCov_Illumina_PE and Augur Workflows](../assets/sops/TG-FLU-PE_SOP_Flu_IlluminaPE_v1.pdf) | v1 | v1.x+ |

### Quality evaluation

The TheiaX workflows will generate various quality metrics. These should be evaluated relative to quality thresholds that have been agreed upon within your laboratory or sequencing program and define the sufficient quality characteristics for a genome and sequence data to be used. For the [TheiaCoV Workflow Series](../workflows/genomic_characterization/theiacov.md), [TheiaProk Workflow Series](../workflows/genomic_characterization/theiaprok.md), and [TheiaEuk Workflow Series](../workflows/genomic_characterization/theiaeuk.md) workflows, this quality evaluation may be undertaken using the optional `QC_check` task. Full instructions for the use of this task may be found on the relevant workflow page. Some quality metrics are not evaluated by the `QC_check` task and should be evaluated manually.

Genomes that fail to meet agreed quality thresholds should not be used. Results for characterization of these genomes may be inaccurate or unreliable. The inclusion of poor-quality genomes in downstream comparative analyses will bias their results. Samples that fail to meet QC thresholds will need to be re-sequenced and sample processing may need to be repeated (e.g. culture-based isolation of clonal bacteria, DNA/RNA extraction, and processing for sequencing).

### Update workflows for SARS-CoV-2 genomes

Workflows are available for updating the Pangolin and VADR assignments made to SARS-CoV-2 genomes. The [Pangolin Update](../workflows/genomic_characterization/pangolin_update.md) workflow accounts for the delay in assigning names to newly emerging lineages that you may have already sequenced. The [VADR_Update](../workflows/genomic_characterization/vadr_update.md) workflow similarly accounts for features that have been newly identified in SARS-CoV-2 genomes when assessing genome quality with VADR.

## Phylogenetics

### Phylogenetic construction

Phylogenetic trees are constructed to assess the evolutionary relationships between sequences in the tree. These evolutionary relationships are often used as a proxy for epidemiological relationships, and sometimes for inferring transmission between isolation sources. 

There are various methods for constructing phylogenetic trees, depending on the sequencing data being used, the organism being analyzed and how it evolved, what you would like to infer from the tree, and the computational resources available for the tree construction. Theiagen has a number of workflows for constructing phylogenetic trees. For full details of these workflows, please see [Guide to Phylogenetics](../guides/phylogenetics.md) which includes advice on the appropriate tree-building workflows and phylogenetic visualization approaches.

!!! example "SOPs for ==phylogenetic construction=="
    | SOP | SOP Version | PHB Version Compatibility |
    |---|---|---|
    | [Analyzing Flu Data in Terra using TheiaCov_Illumina_PE and Augur Workflows](../assets/sops/TG-FLU-PE_SOP_Flu_IlluminaPE_v1.pdf) | v1 | v1.x+ |
    | [Analyzing Phylogenetic Relationships in Terra using Theiagen’s Augur Workflows](../assets/sops/TG-AUGUR-01_Augur_v1.pdf) | v1 | v1.x+ |

### Phylogenetic placement

Phylogenetic placement is used to place your own sequences onto an existing phylogenetic tree. This may be used to find the closest relatives to your sequence(s). More details, including phylogenetic visualization approaches can be found in [Guide to Phylogenetics](../guides/phylogenetics.md).  

## Public Data Sharing

!!! example "SOPs for ==data submissions=="
    | SOP | SOP Version | PHB Version Compatibility |
    |---|---|---|
    | [Submitting SC2 Sequence Data to GISAID using Theiagen’s Terra 2 GISAID Workflow](../assets/sops/TG-GISAID-01_Terra2GISAID_v2.pdf) | v2 | v2+ |

## SARS-CoV-2 Metagenomic Analysis

!!! example "SOPs for ==SARS-CoV-2 metagenomic data analysis=="
    | SOP | SOP Version | PHB Version Compatibility |
    |---|---|---|
    | [Analyzing SARS-CoV-2 Metagenomic Samples using Freyja FASTQ](../assets/sops/TG-FREY-01_FreyjaFASTQ_v2.pdf) | v2 | v2+ |
    | [Plotting SARS-CoV-2 Metagenomic Sample Data using Freyja Plot](../assets/sops/TG-FREY-02_FreyjaPlot_v3.pdf) | v3 | v2+ |
    | [Creating a Dashboard Visualization of SARS-CoV-2 Metagenomic Samples using Freyja Dashboard](../assets/sops/TG-FREY-03-SOP_FreyjaDashboard_v2.pdf) | v2 | v2+ |
    | [Creating Static Reference Files for Freyja Analysis in Terra using Freyja Update](../assets/sops/TG-FREY-04_FreyjaUpdate_v2.pdf) | v2 | v2+ |
