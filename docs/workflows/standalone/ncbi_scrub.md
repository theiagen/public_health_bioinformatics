# NCBI_Scrub

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Standalone](../../workflows_overview/workflows_type.md/#standalone) | [Any Taxa](../../workflows_overview/workflows_kingdom.md/#any-taxa) | PHB v2.2.1 | Yes | Sample-level |

## NCBI Scrub Workflows

NCBI Scrub, also known as the human read removal tool (HRRT), is based on the [SRA Taxonomy Analysis Tool](https://doi.org/10.1186/s13059-021-02490-0) that will take as input a FASTQ file, and produce as output a FASTQ file in which all reads identified as potentially of human origin are either removed (default) or masked with 'N'.
There are three Kraken2 workflows:

- `NCBI_Scrub_PE` is compatible with **Illumina paired-end data**
- `NCBI_Scrub_SE` is compatible with **Illumina single-end data**

### Inputs

<div class="searchable-table" markdown="1">

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** | **Workflow** |
|---|---|---|---|---|---|---|
| dehost_pe or dehost_se | **read1** | File | | | Required | PE, SE |
| dehost_pe or dehost_se | **read2** | File | | | Required | PE |
| dehost_pe or dehost_se | **samplename** | String | | | Required | PE, SE |
| kraken2 | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional | PE, SE |
| kraken2 | **disk_size** | Int | Amount of storage (in GB) to allocate to the task. Increase this when using large (>30GB kraken2 databases such as the "k2_standard" database) | 100 | Optional | PE, SE |
| kraken2 | **docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/kraken2:2.0.8-beta_hv | Optional | PE, SE |
| kraken2 | **kraken2_db** | String | The database used to run Kraken2 | /kraken2-db | Optional | PE, SE |
| kraken2 | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional | PE, SE |
| kraken2 | **read2** | File | Internal component, do not modify | | Do not modify, Optional | SE |
| kraken2 | **target_organism** | String | The organism whose abundance the user wants to check in their reads. This should be a proper taxonomic name recognized by the Kraken database. | | Optional | PE, SE |
| ncbi_scrub_pe or ncbi_scrub_se | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional | PE, SE |
| ncbi_scrub_pe or ncbi_scrub_se | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional | PE, SE |
| ncbi_scrub_pe or  | **docker** | Int | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/ncbi/sra-human-scrubber:2.2.1 | Optional | PE, SE |
| ncbi_scrub_pe or ncbi_scrub_se | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional | PE, SE |
| version_capture | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional | PE, SE |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional | PE, SE |

</div>

### Workflow Tasks

This workflow is composed of two tasks, one to dehost the input reads and another to screen the clean reads with kraken2 and the viral+human database.

??? task "`ncbi_scrub`: human read removal tool"
    Briefly, the HRRT employs a k-mer database constructed of k-mers from Eukaryota derived from all human RefSeq records and subtracts any k-mers found in non-Eukaryota RefSeq records. The remaining set of k-mers compose the database used to identify human reads by the removal tool.

    !!! techdetails "Tool Name Technical Details"
        |  | Links | 
        | --- | --- | 
        | Task | [task_ncbi_scrub.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_ncbi_scrub.wdl) |
        | Software Source Code | [HRRT on GitHub](https://github.com/ncbi/sra-human-scrubber) |
        | Software Documentation | [HRRT on NCBI](https://ncbiinsights.ncbi.nlm.nih.gov/2023/02/02/scrubbing-human-sequences-sra-submissions/) |

??? task "`kraken2`: taxonomic profiling"

    Kraken2 is a bioinformatics tool originally designed for metagenomic applications. It has additionally proven valuable for validating taxonomic assignments and checking contamination of single-species (e.g. bacterial isolate, eukaryotic isolate, viral isolate, etc.) whole genome sequence data.

    Kraken2 is run on the set of raw reads, provided as input, as well as the set of clean reads that are resulted from the `read_QC_trim` workflow

    !!! info "Database-dependent"
        TheiaCoV automatically uses a viral-specific Kraken2 database.

    !!! techdetails "Kraken2 Technical Details"    
        
        |  | Links |
        | --- | --- |
        | Task | [task_kraken2.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/contamination/task_kraken2.wdl) |
        | Software Source Code | [Kraken2 on GitHub](https://github.com/DerrickWood/kraken2/) |
        | Software Documentation | <https://github.com/DerrickWood/kraken2/wiki> |
        | Original Publication(s) | [Improved metagenomic analysis with Kraken 2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0) |

### Outputs

<div class="searchable-table" markdown="1">

| **Variable** | **Type** | **Description** | **Workflow** |
|---|---|---|---|
| kraken_human_dehosted | Float | Percent of human read data detected using the Kraken2 software after host removal | PE, SE |
| kraken_report_dehosted | File | Full Kraken report after host removal | PE, SE |
| kraken_sc2_dehosted | Float | Percent of SARS-CoV-2 read data detected using the Kraken2 software after host removal | PE, SE |
| kraken_version_dehosted | String | Version of Kraken2 software used | PE, SE |
| ncbi_scrub_docker | String | Docker image used to run HRRT | PE, SE |
| ncbi_scrub_human_spots_removed | Int | Number of spots removed (or masked) | PE, SE |
| ncbi_scrub_pe_analysis_date | String | Date of analysis | PE, SE |
| ncbi_scrub_pe_version | String | Version of HRRT software used | PE, SE |
| read1_dehosted | File | Dehosted forward reads | PE, SE |
| read2_dehosted | File | Dehosted reverse reads | PE |

</div>