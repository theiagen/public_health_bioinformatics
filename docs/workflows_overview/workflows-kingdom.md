---
title: Workflows by Kingdom
---

[Sort by Type](workflows-type.md) | [Sort Alphabetically](workflows-alphabetically.md)

---

### Any Taxa


| **Name** | **Description** | **Taxa** | **Workflow type** | **Command-line compatible** | **Last known changes** | **Dockstore** |
|---|---|---|---|---|---|---|
| [**Assembly_Fetch**](../workflows/assembly_fetch.md) | Download assemblies from NCBI, after optionally identifying the closest RefSeq reference genome to your own draft assembly | Any taxa | Sample-level | Yes | v1.3.0 | [Assembly_Fetch_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Assembly_Fetch_PHB:main?tab=info) |
| **BaseSpace_Fetch** | Import data from BaseSpace into Terra | Any taxa | Sample-level | Yes | v2.0.0 | [BaseSpace_Fetch_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/BaseSpace_Fetch_PHB:main?tab=info) |
| **Concatenate_Column_Content** | Concatenate contents of a specified Terra data table column for many samples ("entities") | Any taxa | Set-level | Yes | v2.1.0 | [Concatenate_Column_Content_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Concatenate_Column_Content_PHB:main?tab=info) |
| **Create_Terra_Table** | Upload data to Terra and then run this workflow to have the table automatically created | Any taxa | | Yes | v2.2.0 | [Create_Terra_Table_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Create_Terra_Table_PHB:main?tab=info) |
| **Kraken2** | Taxa identification from reads | Any taxa | Sample-level | Yes | v2.0.0 | [Kraken2_PE_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Kraken2_PE_PHB:main?tab=info), [Kraken2_SE_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Kraken2_SE_PHB:main?tab=info) |
| **RASUSA** | Randomly subsample sequencing reads to a specified coverage | Any taxa | Sample-level | Yes | v2.0.0 | [RASUSA_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/RASUSA_PHB:main?tab=info) |
| **Rename_FASTQ** | Rename paired-end or single-end read files in a Terra data table in a non-destructive way | Any taxa | Sample-level | Yes | v2.1.0 | [Rename_FASTQ_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Rename_FASTQ_PHB:im-utilities-rename-files?tab=info) |
| **SRA_Fetch** | Import publicly available reads from SRA using SRR#, ERR# or DRR# | Any taxa | Sample-level | Yes | v2.2.0 | [SRA_Fetch_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/SRA_Fetch_PHB:main?tab=info) |
| **TheiaMeta** | Genome assembly and QC from metagenomic sequencing | Any taxa | Sample-level | Yes | v2.0.0 | [TheiaMeta_Illumina_PE_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaMeta_Illumina_PE_PHB:main?tab=info) |
| **TheiaValidate** | This workflow performs basic comparisons between user-designated columns in two separate tables. | Any taxa | | No | v2.0.0 | [TheiaValidate_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaValidate_PHB:main?tab=info) |
| **Transfer_Column_Content** | Transfer contents of a specified Terra data table column for many samples (”entities”) to a GCP storage bucket location | Any taxa | Set-level | Yes | v1.3.0 | [Transfer_Column_Content_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Transfer_Column_Content_PHB:main?tab=info) |
| **Zip_Column_Content** | Zip contents of a specified Terra data table column for many samples ("entities") | Any taxa | Set-level | Yes | v2.1.0 | [Zip_Column_Content_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Zip_Column_Content_PHB:main?tab=info) |

### Bacteria

| **Name** | **Description** | **Taxa** | **Workflow type** | **Command-line compatible** | **Last known changes** | **Dockstore** |
|---|---|---|---|---|---|---|
| **Core_Gene_SNP** | Pangenome analysis | Bacteria | Set-level | Some optional features incompatible, Yes | v2.1.0 | [Core_Gene_SNP_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Core_Gene_SNP_PHB:main?tab=info) |
| **Find_Shared_Variants** | Combines and reshapes variant data from Snippy_Variants to illustrate variants shared across multiple samples | Bacteria, Mycotics | Set-level | Yes | v2.0.0 | [Find_Shared_Variants_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Find_Shared_Variants_PHB:main?tab=info) |
| **GAMBIT_Query** | Taxon identification of genome assembly using GAMBIT | Bacteria, Mycotics | Sample-level | Yes | v2.0.0 | [Gambit_Query_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Gambit_Query_PHB:main?tab=info) |
| **kSNP3** | SNP-based phylogenetic analysis from assemblies | Bacteria, Mycotics, Viral | Set-level | Some optional features incompatible, Yes | v2.1.0 | [kSNP3_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/kSNP3_PHB:main?tab=info) |
| **Lyve_SET** | Alignment of reads to a reference genome, SNP calling, curation of high quality SNPs, phylogenetic analysis | Bacteria | Set-level | Yes | v2.1.0 | [Lyve_SET_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Lyve_SET_PHB:main?tab=info) |
| **MashTree_FASTA** | Mash-distance based phylogenetic analysis from assemblies | Bacteria, Mycotics, Viral | Set-level | Some optional features incompatible, Yes | v2.1.0 | [MashTree_FASTA_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/MashTree_FASTA_PHB:main?tab=info) |
| **NCBI-AMRFinderPlus** | Runs NCBI’s AMRFinderPlus on genome assemblies (bacterial and fungal) | Bacteria, Mycotics | Sample-level | Yes | v2.0.0 | [NCBI-AMRFinderPlus_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/NCBI-AMRFinderPlus_PHB:main?tab=info) |
| **Snippy_Streamline** | Implementation of Snippy workflows for phylogenetic analysis from reads, with optional dynamic reference selection | Bacteria | Set-level | Yes | v2.2.0 | [Snippy_Streamline_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Snippy_Streamline_PHB:main?tab=info) |
| **Snippy_Streamline_FASTA** | Implementation of Snippy workflows for phylogenetic analysis from assembled genomes (in FASTA format), with optional dynamic reference selection | Bacteria | Set-level | Yes | v2.2.0 | [Snippy_Streamline_FASTA_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Snippy_Streamline_FASTA_PHB:im-snippy-fasta-dev?tab=info) |
| **Snippy_Tree** | SNP-based phylogenetic analysis from reads, with option to mask recombination | Bacteria | Set-level | Some optional features incompatible, Yes | v2.1.0 | [Snippy_Tree_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Snippy_Tree_PHB:main?tab=info) |
| **Snippy_Variants** | Alignment of reads to a reference genome, then SNP calling | Bacteria, Mycotics, Viral | Sample-level | Yes | v2.2.0 | [Snippy_Variants_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Snippy_Variants_PHB:main?tab=info) |
| **TBProfiler_tNGS** | Performs in silico antimicrobial susceptibility testing on Mycobacterium tuberculosis targeted-NGS samples with TBProfiler and tbp-parser | Bacteria, TB | Sample-level | Yes | v2.0.0 | [TBProfiler_tNGS_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TBProfiler_tNGS_PHB:smw-tngs-tbprofiler-dev?tab=info) |
| **Terra_2_NCBI** | Upload of sequence data to NCBI | Bacteria, Mycotics, Viral | Set-level | No | v2.1.0 | [Terra_2_NCBI_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Terra_2_NCBI_PHB:main?tab=info) |
| **TheiaProk Workflow Series** | Bacterial genome assembly, QC and characterization from WGS data | Bacteria | Sample-level | Some optional features incompatible, Yes | v2.2.0 | [TheiaProk_Illumina_PE_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaProk_Illumina_PE_PHB:main?tab=info), [TheiaProk_Illumina_SE_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaProk_Illumina_SE_PHB:main?tab=info), [TheiaProk_ONT_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaProk_ONT_PHB:main?tab=info), [TheiaProk_FASTA_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaProk_FASTA_PHB:main?tab=info) |

### Mycotics

| **Name** | **Description** | **Taxa** | **Workflow type** | **Command-line compatible** | **Last known changes** | **Dockstore** |
|---|---|---|---|---|---|---|
| **Cauris_CladeTyper** | C. auris clade assignment | Mycotics | Sample-level | Yes | v1.0.0 | [Cauris_CladeTyper_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Cauris_CladeTyper_PHB:main?tab=info) |
| **Find_Shared_Variants** | Combines and reshapes variant data from Snippy_Variants to illustrate variants shared across multiple samples | Bacteria, Mycotics | Set-level | Yes | v2.0.0 | [Find_Shared_Variants_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Find_Shared_Variants_PHB:main?tab=info) |
| **GAMBIT_Query** | Taxon identification of genome assembly using GAMBIT | Bacteria, Mycotics | Sample-level | Yes | v2.0.0 | [Gambit_Query_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Gambit_Query_PHB:main?tab=info) |
| **kSNP3** | SNP-based phylogenetic analysis from assemblies | Bacteria, Mycotics, Viral | Set-level | Some optional features incompatible, Yes | v2.1.0 | [kSNP3_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/kSNP3_PHB:main?tab=info) |
| **MashTree_FASTA** | Mash-distance based phylogenetic analysis from assemblies | Bacteria, Mycotics, Viral | Set-level | Some optional features incompatible, Yes | v2.1.0 | [MashTree_FASTA_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/MashTree_FASTA_PHB:main?tab=info) |
| **NCBI-AMRFinderPlus** | Runs NCBI’s AMRFinderPlus on genome assemblies (bacterial and fungal) | Bacteria, Mycotics | Sample-level | Yes | v2.0.0 | [NCBI-AMRFinderPlus_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/NCBI-AMRFinderPlus_PHB:main?tab=info) |
| **Snippy_Variants** | Alignment of reads to a reference genome, then SNP calling | Bacteria, Mycotics, Viral | Sample-level | Yes | v2.2.0 | [Snippy_Variants_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Snippy_Variants_PHB:main?tab=info) |
| **Terra_2_NCBI** | Upload of sequence data to NCBI | Bacteria, Mycotics, Viral | Set-level | No | v2.1.0 | [Terra_2_NCBI_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Terra_2_NCBI_PHB:main?tab=info) |
| **TheiaEuk** | Mycotic genome assembly, QC and characterization from WGS data | Mycotics | Sample-level | Some optional features incompatible, Yes | v2.0.1 | [TheiaEuk_Illumina_PE_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaEuk_Illumina_PE_PHB:main?tab=info) |

### Viral


| **Name** | **Description** | **Taxa** | **Workflow type** | **Command-line compatible** | **Last known changes** | **Dockstore** |
|---|---|---|---|---|---|---|
| **Augur** | Phylogenetic analysis for viral pathogens | Viral | Sample-level, Set-level | Yes | v2.1.0 | [Augur_Prep_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Augur_Prep_PHB:main?tab=info), [Augur_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Augur_PHB:main?tab=info) |
| **CZGenEpi_Prep** | Prepare metadata and fasta files for easy upload to the CZ GEN EPI platform. | Monkeypox virus, SARS-CoV-2, Viral | Set-level | No | v1.3.0 | [CZGenEpi_Prep_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/CZGenEpi_Prep_PHB:main?tab=info) |
| **Freyja Workflow Series** | Recovers relative lineage abundances from mixed sample data and generates visualizations | SARS-CoV-2, Viral | Sample-level, Set-level | Yes | v2.2.0 | [Freyja_FASTQ_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Freyja_FASTQ_PHB:main?tab=info), [Freyja_Plot_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Freyja_Plot_PHB:main?tab=info), [Freyja_Dashboard_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Freyja_Dashboard_PHB:main?tab=info), [Freyja_Update_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Freyja_Update_PHB:main?tab=info) |
| **kSNP3** | SNP-based phylogenetic analysis from assemblies | Bacteria, Mycotics, Viral | Set-level | Some optional features incompatible, Yes | v2.1.0 | [kSNP3_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/kSNP3_PHB:main?tab=info) |
| **MashTree_FASTA** | Mash-distance based phylogenetic analysis from assemblies | Bacteria, Mycotics, Viral | Set-level | Some optional features incompatible, Yes | v2.1.0 | [MashTree_FASTA_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/MashTree_FASTA_PHB:main?tab=info) |
| **Mercury_Prep_N_Batch** | Prepare metadata and sequence data for submission to NCBI and GISAID | Influenza, Monkeypox virus, SARS-CoV-2, Viral | Set-level | No | v2.2.0 | [Mercury_Prep_N_Batch_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Mercury_Prep_N_Batch_PHB:main?tab=info) |
| **Pangolin Update** | Update Pangolin assignments | SARS-CoV-2, Viral | Sample-level | Yes | v2.0.0 | [Pangolin_Update_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Pangolin_Update_PHB:main?tab=info) |
| **Samples_to_Ref_Tree** | Use Nextclade to rapidly and accurately place your samples on any existing phylogenetic tree | Monkeypox virus, SARS-CoV-2, Viral | Sample-level, Set-level | Yes | v2.1.0 | [Samples_to_Ref_Tree_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Samples_to_Ref_Tree_PHB:main?tab=info) |
| **Snippy_Variants** | Alignment of reads to a reference genome, then SNP calling | Bacteria, Mycotics, Viral | Sample-level | Yes | v2.2.0 | [Snippy_Variants_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Snippy_Variants_PHB:main?tab=info) |
| **Terra_2_GISAID** | Upload of assembly data to GISAID | SARS-CoV-2, Viral | Set-level | Yes | v1.2.1 | [Terra_2_GISAID_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Terra_2_GISAID_PHB:main?tab=info) |
| **Terra_2_NCBI** | Upload of sequence data to NCBI | Bacteria, Mycotics, Viral | Set-level | No | v2.1.0 | [Terra_2_NCBI_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Terra_2_NCBI_PHB:main?tab=info) |
| **TheiaCoV Workflow Series** | Viral genome assembly, QC and characterization from amplicon sequencing | HIV, Influenza, Monkeypox virus, RSV-A, RSV-B, SARS-CoV-2, Viral, WNV | Sample-level, Set-level | Some optional features incompatible, Yes | v2.2.0 | [TheiaCoV_Illumina_PE_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaCoV_Illumina_PE_PHB:main?tab=info), [TheiaCoV_Illumina_SE_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaCoV_Illumina_SE_PHB:main?tab=info), [TheiaCoV_ONT_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaCoV_ONT_PHB:main?tab=info), [TheiaCoV_ClearLabs_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaCoV_ClearLabs_PHB:main?tab=info), [TheiaCoV_FASTA_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaCoV_FASTA_PHB:main?tab=info), [TheiaCoV_FASTA_Batch_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaCoV_FASTA_Batch_PHB:main?tab=info) |
| **Usher_PHB** | Use UShER to rapidly and accurately place your samples on any existing phylogenetic tree | Monkeypox virus, SARS-CoV-2, Viral | Sample-level, Set-level | Yes | v2.1.0 | [Usher_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Usher_PHB:main?tab=info) |
| **VADR_Update** | Update VADR assignments | HAV, Influenza, Monkeypox virus, RSV-A, RSV-B, SARS-CoV-2, Viral, WNV | Sample-level | Yes | v1.2.1 | [VADR_Update_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/VADR_Update_PHB:main?tab=info) |

<!-- definitions for workflow type column -->
*[Sample-level]: This workflow is run once for each sample
*[Set-level]: This workflow is run once on a group of samples

<!-- definitions for taxa column -->
*[Any taxa]: This workflow is organism-agnostic and can be run with any taxa
*[Viral]: This workflow is compatible with any viral pathogen
*[Bacteria]: This workflow is compatible with any bacterial pathogen
*[Mycotics]: This workflow is compatible with mycotic pathogens

<!-- definition for command-line compatible column -->
*[Command-line compatible]: Command-line compatibility is determined if the workflow can be run on a local command-line environment, providing all dependencies are installed, with either `miniwdl` or `cromwell`.
*[Some optional features incompatible]: Some optional features of this workflow are incompatible with command-line use and require modification
*[Yes]: This workflow is compatible with command-line use
*[No]: This workflow is not compatible with command-line use even with modifications
