# AMR Search

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filter_column="Name", filter_values="[**AMR_Search**](../workflows/standalone/amr_search.md)", columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level"]) }}

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

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="AMR_Search", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

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

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filter_column="Workflow", filter_values="AMR_Search", columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///

## References

> [Pathogenwatch AMRsearch](https://github.com/pathogenwatch-oss/amr-search)
<!-- -->
> [PAARSNP](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7545138/)