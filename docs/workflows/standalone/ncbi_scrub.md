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

!!! caption ""
    === "NCBI_Scrub_PE"
        /// html | div[class="searchable-table"]

        {{ input_table("docs/assets/input_tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="NCBI_Scrub_PE", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"], indent=8) }}
        ///

    === "NCBI_Scrub_SE"
        /// html | div[class="searchable-table"]

        {{ input_table("docs/assets/input_tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="NCBI_Scrub_SE", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"], indent=8) }}
        ///

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
