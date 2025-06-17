# AMR Search

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Standalone](../../workflows_overview/workflows_type.md/#standalone) | [Any Taxa](../../workflows_overview/workflows_kingdom.md/#any-taxa) | v3.0.0 | Yes | Sample-level |

## AMR_Search_PHB

!!! caption "AMR_Search Workflow Diagram"
    ![AMR_Search Workflow Diagram](../../assets/figures/AMR_Search.png)

The AMR_Search workflow is a standalone version of Pathogenwatch's AMR profiling functionality utilizing `AMRsearch` tool from Pathogenwatch.

A limited number of species are currently supported and are listed below. NCBI codes are needed from this table to select the correct library.

| Species                      | NCBI Code |
|------------------------------|-----------|
| _Neisseria gonorrhoeae_      | 485       |
| _Staphylococcus aureus_      | 1280      |
| _Salmonella typhi_          | 90370     |
| _Streptococcus pneumoniae_   | 1313      |
| _Klebisiella_                | 570       |
| _Klebsiella pneumoniae_     | 573       |
| _Candida auris_              | 498019    |
| _Vibrio cholerae_            | 666       |
| _Campylobacter_              | 194       |

### Inputs

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| amr_search_workflow | **amr_search_database** | String | NCBI taxon code of samples known taxonomy, see above supported species || Required |
| amr_search_workflow | **input_fasta** | File | A microbial assembly file || Required |
| amr_search_workflow | **samplename** | String | Identifier user wants prefixed to output files || Required |
| amr_search | **cpu** | Integer | Number of CPUs to allocate to the task |2| Optional |
| amr_search | **disk_size** | Integer | Amount of storage (in GB) to allocate to the task |50| Optional |
| amr_search | **docker** | String | The docker container to use for the task |us-docker.pkg.dev/general-theiagen/theiagen/amrsearch:0.2.0| Optional |
| amr_search | **memory** | Integer | Amount of memory/RAM (in GB) to allocate to the task |8| Optional |

### Workflow Tasks

??? task "`amr_search`: Antimicrobial resistance profiling"

    This task performs *in silico* antimicrobial resistance (AMR) profiling for supported species using **AMRsearch**, the primary tool used by [Pathogenwatch](https://pathogen.watch/) to genotype and infer antimicrobial resistance (AMR) phenotypes from assembled microbial genomes.

    **AMRsearch** screens against Pathogenwatch's library of curated genotypes and inferred phenotypes, developed in collaboration with community experts. Resistance phenotypes are determined based on both **resistance genes** and **mutations**, and the system accounts for interactions between multiple SNPs, genes, and suppressors. Predictions follow **S/I/R classification** (*Sensitive, Intermediate, Resistant*).

    **Outputs:**

    - **JSON Output:** Contains the complete AMR profile, including detailed **resistance state**, detected **resistance genes/mutations**, and supporting **BLAST results**.

    - **CSV & PDF Tables:** An incorprated Python script, `parse_amr_json.py`, extracts and formats results into a **CSV file** and **PDF summary table** for easier visualization.

    !!! techdetails "amr_search Technical Details"    

        |  | Links |
        | --- | --- |
        | Task | [task_amr_search.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/drug_resistance/task_amr_search.wdl) |
        | Software Source Code | [AMRsearch](https://github.com/pathogenwatch-oss/amr-search) |
        | Software Documentation | [Pathogenwatch](https://cgps.gitbook.io/pathogenwatch) |
        | Original Publication(s) | [PAARSNP: *rapid genotypic resistance prediction for *Neisseria gonorrhoeae*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7545138/) |

### Outputs

| **Variable** | **Type** | **Description** |
|---|---|---|
| amr_results_csv | File | CSV formatted AMR profile |
| amr_results_pdf | File | PDF formatted AMR profile |
| amr_search_results | File | JSON formatted AMR profile including BLAST results |
| amr_search_docker | String | Docker image used to run AMR_Search |
| amr_search_version | String | Version of AMR_Search libraries used |

## References

> [Pathogenwatch AMRsearch](https://github.com/pathogenwatch-oss/amr-search)
<!-- -->
> [PAARSNP](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7545138/)