---
title: Workflows by Type
---

[Sort by Kingdom](workflows_kingdom.md) | [Sort Alphabetically](workflows_alphabetically.md)

---

### Data Import

<div class="searchable-table" markdown="1">

| **Name** | **Description** | **Applicable Kingdom** | **Workflow Level** | **Command-line Compatibility**[^1] | **Last Known Changes** | **Dockstore** |
|---|---|---|---|---|---|---|
| [**Assembly_Fetch**](../workflows/data_import/assembly_fetch.md) | Download assemblies from NCBI, after optionally identifying the closest RefSeq reference genome to your own draft assembly | Any taxa | Sample-level | Yes | v1.3.0 | [Assembly_Fetch_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Assembly_Fetch_PHB:main?tab=info) |
| [**BaseSpace_Fetch**](../workflows/data_import/basespace_fetch.md)| Import data from BaseSpace into Terra | Any taxa | Sample-level | Yes | v2.0.0 | [BaseSpace_Fetch_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/BaseSpace_Fetch_PHB:main?tab=info) |
| [**Create_Terra_Table**](../workflows/data_import/create_terra_table.md)| Upload data to Terra and then run this workflow to have the table automatically created | Any taxa | | Yes | v2.2.0 | [Create_Terra_Table_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Create_Terra_Table_PHB:main?tab=info) |
| [**SRA_Fetch**](../workflows/data_import/sra_fetch.md)| Import publicly available reads from SRA using SRR#, ERR# or DRR# | Any taxa | Sample-level | Yes | v2.2.0 | [SRA_Fetch_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/SRA_Fetch_PHB:main?tab=info) |

</div>

### Genomic Characterization

<div class="searchable-table" markdown="1">

| **Name** | **Description** | **Applicable Kingdom** | **Workflow Level** | **Command-line Compatibility**[^1] | **Last Known Changes** | **Dockstore** |
|---|---|---|---|---|---|---|
| [**Freyja Workflow Series**](../workflows/genomic_characterization/freyja.md)| Recovers relative lineage abundances from mixed sample data and generates visualizations | SARS-CoV-2, Viral | Sample-level, Set-level | Yes | v2.2.0 | [Freyja_FASTQ_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Freyja_FASTQ_PHB:main?tab=info), [Freyja_Plot_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Freyja_Plot_PHB:main?tab=info), [Freyja_Dashboard_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Freyja_Dashboard_PHB:main?tab=info), [Freyja_Update_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Freyja_Update_PHB:main?tab=info) |
| [**Pangolin_Update**](../workflows/genomic_characterization/pangolin_update.md) | Update Pangolin assignments | SARS-CoV-2, Viral | Sample-level | Yes | v2.0.0 | [Pangolin_Update_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Pangolin_Update_PHB:main?tab=info) |
| [**TheiaCov Workflow Series**](../workflows/genomic_characterization/theiacov.md) | Viral genome assembly, QC and characterization from amplicon sequencing | HIV, Influenza, Monkeypox virus, RSV-A, RSV-B, SARS-CoV-2, Viral, WNV | Sample-level, Set-level | Some optional features incompatible, Yes | v2.2.0 | [TheiaCoV_Illumina_PE_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaCoV_Illumina_PE_PHB:main?tab=info), [TheiaCoV_Illumina_SE_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaCoV_Illumina_SE_PHB:main?tab=info), [TheiaCoV_ONT_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaCoV_ONT_PHB:main?tab=info), [TheiaCoV_ClearLabs_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaCoV_ClearLabs_PHB:main?tab=info), [TheiaCoV_FASTA_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaCoV_FASTA_PHB:main?tab=info), [TheiaCoV_FASTA_Batch_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaCoV_FASTA_Batch_PHB:main?tab=info) |
| [**TheiaEuk**](../workflows/genomic_characterization/theiaeuk.md) | Mycotic genome assembly, QC and characterization from WGS data | Mycotics | Sample-level | Some optional features incompatible, Yes | v2.0.1 | [TheiaEuk_Illumina_PE_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaEuk_Illumina_PE_PHB:main?tab=info) |
| [**TheiaMeta**](../workflows/genomic_characterization/theiameta.md) | Genome assembly and QC from metagenomic sequencing | Any taxa | Sample-level | Yes | v2.0.0 | [TheiaMeta_Illumina_PE_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaMeta_Illumina_PE_PHB:main?tab=info) |
| [**TheiaProk Workflow Series**](../workflows/genomic_characterization/theiaprok.md) | Bacterial genome assembly, QC and characterization from WGS data | Bacteria | Sample-level | Some optional features incompatible, Yes | v2.2.0 | [TheiaProk_Illumina_PE_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaProk_Illumina_PE_PHB:main?tab=info), [TheiaProk_Illumina_SE_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaProk_Illumina_SE_PHB:main?tab=info), [TheiaProk_ONT_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaProk_ONT_PHB:main?tab=info), [TheiaProk_FASTA_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaProk_FASTA_PHB:main?tab=info) |
| [**VADR_Update**](../workflows/genomic_characterization/vadr_update.md)| Update VADR assignments | HAV, Influenza, Monkeypox virus, RSV-A, RSV-B, SARS-CoV-2, Viral, WNV | Sample-level | Yes | v1.2.1 | [VADR_Update_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/VADR_Update_PHB:main?tab=info) |

</div>

### Phylogenetic Construction

<div class="searchable-table" markdown="1">

| **Name** | **Description** | **Applicable Kingdom** | **Workflow Level** | **Command-line Compatibility**[^1] | **Last Known Changes** | **Dockstore** |
|---|---|---|---|---|---|---|
| [**Augur**](../workflows/phylogenetic_construction/augur.md) | Phylogenetic analysis for viral pathogens | Viral | Sample-level, Set-level | Yes | v2.1.0 | [Augur_Prep_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Augur_Prep_PHB:main?tab=info), [Augur_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Augur_PHB:main?tab=info) |
| [**Core_Gene_SNP**](../workflows/phylogenetic_construction/core_gene_snp.md) | Pangenome analysis | Bacteria | Set-level | Some optional features incompatible, Yes | v2.1.0 | [Core_Gene_SNP_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Core_Gene_SNP_PHB:main?tab=info) |
| [**CZGenEpi_Prep**](../workflows/phylogenetic_construction/czgenepi_prep.md)| Prepare metadata and fasta files for easy upload to the CZ GEN EPI platform. | Monkeypox virus, SARS-CoV-2, Viral | Set-level | No | v1.3.0 | [CZGenEpi_Prep_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/CZGenEpi_Prep_PHB:main?tab=info) |
| [**Find_Shared_Variants**](../workflows/phylogenetic_construction/find_shared_variants.md)| Combines and reshapes variant data from Snippy_Variants to illustrate variants shared across multiple samples | Bacteria, Mycotics | Set-level | Yes | v2.0.0 | [Find_Shared_Variants_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Find_Shared_Variants_PHB:main?tab=info) |
| [**kSNP3**](../workflows/phylogenetic_construction/ksnp3.md)| SNP-based phylogenetic analysis from assemblies | Bacteria, Mycotics, Viral | Set-level | Some optional features incompatible, Yes | v2.1.0 | [kSNP3_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/kSNP3_PHB:main?tab=info) |
| [**Lyve_SET**](../workflows/phylogenetic_construction/lyve_set.md)| Alignment of reads to a reference genome, SNP calling, curation of high quality SNPs, phylogenetic analysis | Bacteria | Set-level | Yes | v2.1.0 | [Lyve_SET_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Lyve_SET_PHB:main?tab=info) |
| [**MashTree_FASTA**](../workflows/phylogenetic_construction/mashtree_fasta.md)| Mash-distance based phylogenetic analysis from assemblies | Bacteria, Mycotics, Viral | Set-level | Some optional features incompatible, Yes | v2.1.0 | [MashTree_FASTA_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/MashTree_FASTA_PHB:main?tab=info) |
| [**Snippy_Streamline**](../workflows/phylogenetic_construction/snippy_streamline.md)| Implementation of Snippy workflows for phylogenetic analysis from reads, with optional dynamic reference selection | Bacteria | Set-level | Yes | v2.2.0 | [Snippy_Streamline_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Snippy_Streamline_PHB:main?tab=info) |
| [**Snippy_Streamline_FASTA**](../workflows/phylogenetic_construction/snippy_streamline_fasta.md)| Implementation of Snippy workflows for phylogenetic analysis from assembled genomes (in FASTA format), with optional dynamic reference selection | Bacteria | Set-level | Yes | v2.2.0 | [Snippy_Streamline_FASTA_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Snippy_Streamline_FASTA_PHB:im-snippy-fasta-dev?tab=info) |
| [**Snippy_Tree**](../workflows/phylogenetic_construction/snippy_tree.md)| SNP-based phylogenetic analysis from reads, with option to mask recombination | Bacteria | Set-level | Some optional features incompatible, Yes | v2.1.0 | [Snippy_Tree_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Snippy_Tree_PHB:main?tab=info) |
| [**Snippy_Variants**](../workflows/phylogenetic_construction/snippy_variants.md)| Alignment of reads to a reference genome, then SNP calling | Bacteria, Mycotics, Viral | Sample-level | Yes | v2.2.0 | [Snippy_Variants_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Snippy_Variants_PHB:main?tab=info) |

</div>

### Phylogenetic Placement

<div class="searchable-table" markdown="1">

| **Name** | **Description** | **Applicable Kingdom** | **Workflow Level** | **Command-line Compatibility**[^1] | **Last Known Changes** | **Dockstore** |
|---|---|---|---|---|---|---|
| [**Samples_to_Ref_Tree**](../workflows/phylogenetic_placement/samples_to_ref_tree.md)| Use Nextclade to rapidly and accurately place your samples on any existing phylogenetic tree | Monkeypox virus, SARS-CoV-2, Viral | Sample-level, Set-level | Yes | v2.1.0 | [Samples_to_Ref_Tree_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Samples_to_Ref_Tree_PHB:main?tab=info) |
| [**Usher_PHB**](../workflows/phylogenetic_placement/usher.md)| Use UShER to rapidly and accurately place your samples on any existing phylogenetic tree | Monkeypox virus, SARS-CoV-2, Viral | Sample-level, Set-level | Yes | v2.1.0 | [Usher_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Usher_PHB:main?tab=info) |

</div>

### Public Data Sharing

<div class="searchable-table" markdown="1">

| **Name** | **Description** | **Applicable Kingdom** | **Workflow Level** | **Command-line Compatibility**[^1] | **Last Known Changes** | **Dockstore** |
|---|---|---|---|---|---|---|
| [**Mercury_Prep_N_Batch**](../workflows/public_data_sharing/mercury_prep_n_batch.md)| Prepare metadata and sequence data for submission to NCBI and GISAID | Influenza, Monkeypox virus, SARS-CoV-2, Viral | Set-level | No | v2.2.0 | [Mercury_Prep_N_Batch_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Mercury_Prep_N_Batch_PHB:main?tab=info) |
| [**Terra_2_GISAID**](../workflows/public_data_sharing/terra_2_gisaid.md)| Upload of assembly data to GISAID | SARS-CoV-2, Viral | Set-level | Yes | v1.2.1 | [Terra_2_GISAID_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Terra_2_GISAID_PHB:main?tab=info) |
| [**Terra_2_NCBI**](../workflows/public_data_sharing/terra_2_ncbi.md)| Upload of sequence data to NCBI | Bacteria, Mycotics, Viral | Set-level | No | v2.1.0 | [Terra_2_NCBI_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Terra_2_NCBI_PHB:main?tab=info) |

</div>

### Exporting Data from Terra

<div class="searchable-table" markdown="1">

| **Name** | **Description** | **Applicable Kingdom** | **Workflow Level** | **Command-line Compatibility**[^1] | **Last Known Changes** | **Dockstore** |
|---|---|---|---|---|---|---|
| [**Concatenate_Column_Content**](../workflows/data_export/concatenate_column_content.md) | Concatenate contents of a specified Terra data table column for many samples ("entities") | Any taxa | Set-level | Yes | v2.1.0 | [Concatenate_Column_Content_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Concatenate_Column_Content_PHB:main?tab=info) |
| [**Transfer_Column_Content**](../workflows/data_export/transfer_column_content.md)| Transfer contents of a specified Terra data table column for many samples ("entities") to a GCP storage bucket location | Any taxa | Set-level | Yes | v1.3.0 | [Transfer_Column_Content_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Transfer_Column_Content_PHB:main?tab=info) |
| [**Zip_Column_Content**](../workflows/data_export/zip_column_content.md)| Zip contents of a specified Terra data table column for many samples ("entities") | Any taxa | Set-level | Yes | v2.1.0 | [Zip_Column_Content_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Zip_Column_Content_PHB:main?tab=info) |

</div>

### Standalone

<div class="searchable-table" markdown="1">

| **Name** | **Description** | **Applicable Kingdom** | **Workflow Level** | **Command-line Compatibility**[^1] | **Last Known Changes** | **Dockstore** |
|---|---|---|---|---|---|---|
| [**Cauris_CladeTyper**](../workflows/standalone/cauris_cladetyper.md)| C. auris clade assignment | Mycotics | Sample-level | Yes | v1.0.0 | [Cauris_CladeTyper_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Cauris_CladeTyper_PHB:main?tab=info) |
| [**GAMBIT_Query**](../workflows/standalone/gambit_query.md)| Taxon identification of genome assembly using GAMBIT | Bacteria, Mycotics | Sample-level | Yes | v2.0.0 | [Gambit_Query_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Gambit_Query_PHB:main?tab=info) |
| [**Kraken2**](../workflows/standalone/kraken2.md) | Taxa identification from reads | Any taxa | Sample-level | Yes | v2.0.0 | [Kraken2_PE_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Kraken2_PE_PHB:main?tab=info), [Kraken2_SE_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Kraken2_SE_PHB:main?tab=info) |
| [**NCBI-AMRFinderPlus**](../workflows/standalone/ncbi_amrfinderplus.md)| Runs NCBI's AMRFinderPlus on genome assemblies (bacterial and fungal) | Bacteria, Mycotics | Sample-level | Yes | v2.0.0 | [NCBI-AMRFinderPlus_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/NCBI-AMRFinderPlus_PHB:main?tab=info) |
| [**NCBI_Scrub**](../workflows/standalone/ncbi_scrub.md)| Runs NCBI's HRRT on Illumina FASTQs | Any taxa | Sample-level | Yes | v2.2.1 | [NCBI_Scrub_PE_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/NCBI_Scrub_PE_PHB:main?tab=info), [NCBI_Scrub_SE_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/NCBI_Scrub_SE_PHB:main?tab=info) |
| [**RASUSA**](../workflows/standalone/rasusa.md)| Randomly subsample sequencing reads to a specified coverage | Any taxa | Sample-level | Yes | v2.0.0 | [RASUSA_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/RASUSA_PHB:main?tab=info) |
| [**Rename_FASTQ**](../workflows/standalone/rename_fastq.md)| Rename paired-end or single-end read files in a Terra data table in a non-destructive way | Any taxa | Sample-level | Yes | v2.1.0 | [Rename_FASTQ_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/Rename_FASTQ_PHB:im-utilities-rename-files?tab=info) |
| [**TBProfiler_tNGS**](../workflows/standalone/tbprofiler_tngs.md)| Performs in silico antimicrobial susceptibility testing on Mycobacterium tuberculosis targeted-NGS samples with TBProfiler and tbp-parser | Bacteria, TB | Sample-level | Yes | v2.0.0 | [TBProfiler_tNGS_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TBProfiler_tNGS_PHB:smw-tngs-tbprofiler-dev?tab=info) |
| [**TheiaValidate**](../workflows/standalone/theiavalidate.md)| This workflow performs basic comparisons between user-designated columns in two separate tables. | Any taxa | | No | v2.0.0 | [TheiaValidate_PHB](https://dockstore.org/workflows/github.com/theiagen/public_health_bioinformatics/TheiaValidate_PHB:main?tab=info) |

</div>

<!-- definitions for workflow type column -->
*[Sample-level]: This workflow is run once for each sample
*[Set-level]: This workflow is run once on a group of samples

<!-- definitions for taxa column -->
*[Any taxa]: This workflow is organism-agnostic and can be run with any taxa
*[Viral]: This workflow is compatible with any viral pathogen
*[Bacteria]: This workflow is compatible with any bacterial pathogen
*[Mycotics]: This workflow is compatible with mycotic pathogens

<!-- definition for command-line compatible column -->
[^1]:
    Command-line compatibility is determined if the workflow can be run on a local command-line environment, providing all dependencies are installed, with either `miniwdl` or `cromwell`.
*[Some optional features incompatible]: Some optional features of this workflow are incompatible with command-line use and require modification
*[Yes]: This workflow is compatible with command-line use
*[No]: This workflow is not compatible with command-line use even with modifications
