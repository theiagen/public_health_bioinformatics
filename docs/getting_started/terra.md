# Getting Started with Genomic Analysis in Public Health

<aside>
<img src="https://prod-files-secure.s3.us-west-2.amazonaws.com/be290196-9090-4f3c-a9ab-fe730ad213e0/b61fcdac-c45a-44d3-86b9-a241232128cc/download.png" alt="https://prod-files-secure.s3.us-west-2.amazonaws.com/be290196-9090-4f3c-a9ab-fe730ad213e0/b61fcdac-c45a-44d3-86b9-a241232128cc/download.png" width="40px" /> **Theiagen’s approach to genomic analysis in public health typically uses the [Terra](https://terra.bio/) platform to run workflows that undertake bioinformatic analysis, then uses other platforms for visualization of the resulting data. This is described in more depth in our paper “[Accelerating bioinformatics implementation in public health](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001051?crawler=true)”, and the application of this approach for genomic surveillance of SARS-CoV-2 in California is described in the paper “[Pathogen genomics in public health laboratories: successes, challenges, and lessons learned from California’s SARS-CoV-2 Whole-Genome Sequencing Initiative, California COVIDNet](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001027?utm_source=TrendMD&utm_medium=cpc&utm_campaign=Microbial_genomics_TrendMD_0#)”.**

</aside>

<aside>
<img src="https://prod-files-secure.s3.us-west-2.amazonaws.com/be290196-9090-4f3c-a9ab-fe730ad213e0/39a844d4-3053-4f0c-a8c3-dc265e8f9325/Picture3.png" alt="https://prod-files-secure.s3.us-west-2.amazonaws.com/be290196-9090-4f3c-a9ab-fe730ad213e0/39a844d4-3053-4f0c-a8c3-dc265e8f9325/Picture3.png" width="40px" /> **When undertaking genomic analysis using Terra and other data visualization platforms, it is essential to consider the necessary and appropriate workflows and resources for your analysis. To help you make these choices, take a look at the relationship between the most commonly used Theiagen workflows** [in the diagram](https://www.notion.so/Guide-to-Getting-Started-with-Genomic-Analysis-in-Public-Health-3a6dde9a07ef471a9c0bcda6edc52d6a?pvs=21)**, and the descriptions of the major stages in genomic data analysis below.**

</aside>

<aside>
<img src="https://prod-files-secure.s3.us-west-2.amazonaws.com/be290196-9090-4f3c-a9ab-fe730ad213e0/a7a39960-0058-472b-bafa-5109dd1bd393/Picture3.png" alt="https://prod-files-secure.s3.us-west-2.amazonaws.com/be290196-9090-4f3c-a9ab-fe730ad213e0/a7a39960-0058-472b-bafa-5109dd1bd393/Picture3.png" width="40px" /> **Detailed documentation for each PHB release, including helpful workflow input and output explanations, can be found on the Public Health Resources page!**

[**Theiagen Public Health Resources**](https://www.notion.so/Theiagen-Public-Health-Resources-a4bd134b0c5c4fe39870e21029a30566?pvs=21)

</aside>

![**Analysis Approaches for Genomic Data:** This diagram shows the Theiagen workflows (green boxes) available on Terra for analysis of genomic data in public health and the workflows that may be used consecutively (arrows). The blue boxes describe the major functions that these workflows undertake. Descriptions of these functions and their workflows can be found in the [Getting Started with Genomic Analysis in Public Health](https://www.notion.so/PHB-main-DO-NOT-EDIT-7ffec467deea4bcca85576040cad6262?pvs=21) section below. The yellow boxes show functions that may be undertaken independently of workflows on Terra.](https://prod-files-secure.s3.us-west-2.amazonaws.com/be290196-9090-4f3c-a9ab-fe730ad213e0/d16a913d-2aa6-46af-a3f7-fd14390b47f0/Workflow__Relationships_(4)_(1).png)

**Analysis Approaches for Genomic Data:** This diagram shows the Theiagen workflows (green boxes) available on Terra for analysis of genomic data in public health and the workflows that may be used consecutively (arrows). The blue boxes describe the major functions that these workflows undertake. Descriptions of these functions and their workflows can be found in the [Getting Started with Genomic Analysis in Public Health](https://www.notion.so/PHB-main-DO-NOT-EDIT-7ffec467deea4bcca85576040cad6262?pvs=21) section below. The yellow boxes show functions that may be undertaken independently of workflows on Terra.

# Data Import to Terra

To start using Terra for data analysis, you will first need to import your data into your workspace. There are multiple ways to do this:

- **Using Terra’s native features to upload data from your local computer or link to data that’s already in a Google bucket**
- Data import workflows
    - Using the [SRA_Fetch](https://www.notion.so/SRA_Fetch-46bef1c62593438b9a6c8b9256442fa1?pvs=21) workflow to import publicly available data from any repository in the [INSDC](https://www.insdc.org/) (including with [SRA](https://www.ncbi.nlm.nih.gov/sra), [ENA](https://www.ebi.ac.uk/ena/browser/home) and [DRA](https://www.ddbj.nig.ac.jp/dra/index-e.html))
    - Using the [Assembly_Fetch ](https://www.notion.so/Assembly_Fetch-fe838a9ada184166975c6ff7c143e99b?pvs=21) workflow to import publicly available genome assemblies from [NCBI](https://www.ncbi.nlm.nih.gov/datasets/)
    - Using the [BaseSpace_Fetch](https://www.notion.so/BaseSpace_Fetch-ebece1219552422cb143cb316ad59b5b?pvs=21) workflow to import data from your [Illumina BaseSpace](https://basespace.illumina.com/) account
    
    <aside>
    ✅ SOPs
    
    - SOPs for **importing data into a Terra workspace**
        
        [Untitled](https://www.notion.so/4027ff128dd34b13b68757365fc45e7f?pvs=21)
        
    </aside>
    

# Genome assembly, QC, and characterization

## Theia workflows

The Theia workflows are used for genome assembly, quality control, and characterization. The [TheiaCoV Workflow Series](https://www.notion.so/TheiaCoV-Workflow-Series-47e3269a932343be982507d72fcb0fbe?pvs=21), [TheiaProk Workflow Series](https://www.notion.so/TheiaProk-Workflow-Series-6ae7afc7b3544b889fda7c758ea436e1?pvs=21), and [TheiaEuk Workflow Series](https://www.notion.so/TheiaEuk-Workflow-Series-516a773fb418424c876792ea8dba3fde?pvs=21) workflows are intended for viral, bacterial, and fungal pathogens, respectively. [TheiaMeta Workflow Series](https://www.notion.so/TheiaMeta-Workflow-Series-a9a9d5d52e02474189c1d7b3225f3ffe?pvs=21)  is intended for the analysis of a single taxon from metagenomic data.

<aside>
✅ **SOPs**

- SOPs for analyzing **SARS-COV-2** using the TheiaCoV workflows
    
    [Untitled](https://www.notion.so/f03553132dd34a7a95ff2fa540554445?pvs=21)
    
- SOPs for analyzing **influenza** using the TheiaCoV and Augur workflows
    
    [Untitled](https://www.notion.so/4c9c3400945548c98a1e60260b6c8597?pvs=21)
    
</aside>

## Quality evaluation

The TheiaX workflows will generate various quality metrics. These should be evaluated relative to quality thresholds that have been agreed upon within your laboratory or sequencing program and define the sufficient quality characteristics for a genome and sequence data to be used. For the [TheiaCoV Workflow Series](https://www.notion.so/TheiaCoV-Workflow-Series-47e3269a932343be982507d72fcb0fbe?pvs=21), [TheiaProk Workflow Series](https://www.notion.so/TheiaProk-Workflow-Series-6ae7afc7b3544b889fda7c758ea436e1?pvs=21), and [TheiaEuk Workflow Series](https://www.notion.so/TheiaEuk-Workflow-Series-516a773fb418424c876792ea8dba3fde?pvs=21) workflows, this quality evaluation may be undertaken using the optional `QC_check` task. Full instructions for the use of this task may be found on the relevant workflow page. Some quality metrics are not evaluated by the `QC_check` task and should be evaluated manually. 

Genomes that fail to meet agreed quality thresholds should not be used. Results for characterization of these genomes may be inaccurate or unreliable. The inclusion of poor-quality genomes in downstream comparative analyses will bias their results. Samples that fail to meet QC thresholds will need to be re-sequenced and sample processing may need to be repeated (e.g. culture-based isolation of clonal bacteria, DNA/RNA extraction, and processing for sequencing).

## Update workflows for SARS-CoV-2 genomes

Workflows are available for updating the Pangolin and VADR assignments made to SARS-CoV-2 genomes. The [Pangolin Update](https://www.notion.so/Pangolin-Update-4763e30105f343a5ab6a1a1dd315a976?pvs=21)  workflow accounts for the delay in assigning names to newly emerging lineages that you may have already sequenced. The [VADR_Update](https://www.notion.so/VADR_Update-a6398bceef7f408b8b3a7334a2b243c5?pvs=21)  workflow similarly accounts for features that have been newly identified in SARS-CoV-2 genomes when assessing genome quality with VADR.

# Phylogenetics

## Phylogenetic construction

Phylogenetic trees are constructed to assess the evolutionary relationships between sequences in the tree. These evolutionary relationships are often used as a proxy for epidemiological relationships, and sometimes for inferring transmission between isolation sources. 

There are various methods for constructing phylogenetic trees, depending on the sequencing data being used, the organism being analyzed and how it evolved, what you would like to infer from the tree, and the computational resources available for the tree construction. Theiagen has a number of workflows for constructing phylogenetic trees. For full details of these workflows, please see [Guide to Phylogenetics](https://www.notion.so/Guide-to-Phylogenetics-c997fe59e3f0423aa8a73eeccccd1b92?pvs=21) which includes advice on the appropriate tree-building workflows and phylogenetic visualization approaches.

<aside>
✅ **SOPs**

- SOPs for **phylogenetic construction** using Theiagen workflows
    
    [Untitled](https://www.notion.so/dc2e1bc49fba4fb8b710046c87f29bd2?pvs=21)
    
</aside>

## Phylogenetic placement

Phylogenetic placement is used to place your own sequences onto an existing phylogenetic tree. This may be used to find the closest relatives to your sequence(s). More details, including phylogenetic visualization approaches can be found in [Guide to Phylogenetics](https://www.notion.so/Guide-to-Phylogenetics-c997fe59e3f0423aa8a73eeccccd1b92?pvs=21)  

# Public Data Sharing

<aside>
✅ SOPs

- SOPs for submission of **SARS-CoV-2** data to **GISAID** using the Terra_2_GISAID workflow
    
    [Untitled](https://www.notion.so/ae71a26243624e84b0d91f1f431cde79?pvs=21)
    
</aside>

# **SARS-CoV-2 Metagenomic Analysis**

<aside>
✅ SOPs

- SOPs for analyzing **SARS-CoV-2 metagenomic data** using the Freyja workflows
    
    [Untitled](https://www.notion.so/8ecc3d32b79d4c5b8f7521d47d6d5703?pvs=21)
    
</aside>