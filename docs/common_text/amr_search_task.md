<!-- if: amr_search -->
??? task "`amr_search`: Antimicrobial resistance profiling"
<!-- endif -->

<!-- if: theiaprok|theiaeuk -->
??? task "`amr_search`: Antimicrobial Resistance Profiling (optional)"
    To activate this task, set `run_amr_search` to be `true`.
<!-- endif -->

    This task performs _in silico_ antimicrobial resistance (AMR) profiling for supported species using AMRsearch, the primary tool used by [Pathogenwatch](https://pathogen.watch/) to genotype and infer antimicrobial resistance (AMR) phenotypes from assembled microbial genomes.

    AMRsearch screens against Pathogenwatch's library of curated genotypes and inferred phenotypes, developed in collaboration with community experts. Resistance phenotypes are determined based on both _resistance genes_ and _mutations_, and the system accounts for interactions between multiple SNPs, genes, and suppressors. Predictions follow **S/I/R classification** (_Sensitive_, _Intermediate_, _Resistant_).

    Currently, only a subset of species are supported by this task.

    ??? toggle "Supported Species"
        The following table shows the species name and the associated NCBI Code. If you are running AMR Search as part of TheiaProk and TheiaEuk, these codes will be automatically determined based on the GAMBIT predicted taxon, or the user-provided `expected_taxon` input.

        | Species                      | NCBI Code |
        |------------------------------|-----------|
        | _Neisseria gonorrhoeae_      | 485       |
        | _Staphylococcus aureus_      | 1280      |
        | _Salmonella_ Typhi           | 90370     |
        | _Streptococcus pneumoniae_   | 1313      |
        | _Klebisiella_                | 570       |
        | _Escherichia_                | 561       |
        | _Mycobacterium tuberculosis_ | 1773      |
        | _Candida auris_              | 498019    |
        | _Vibrio cholerae_            | 666       |
        | _Campylobacter_              | 194       |

    **Outputs:**

    - **JSON Output**: Contains the complete AMR profile, including detailed resistance state, detected resistance genes/mutations, and supporting BLAST results.
    - **CSV & PDF Tables**: An incorporated Python script, `parse_amr_json.py`, extracts and formats results into a CSV file and PDF summary table for easier visualization.

    !!! techdetails "amr_search Technical Details"    
        |  | Links |
        | --- | --- |
        | Task | [task_amr_search.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/drug_resistance/task_amr_search.wdl) |
        | Software Source Code | [AMRsearch on GitHub](https://github.com/pathogenwatch-oss/amr-search) |
        | Software Documentation | [AMRsearch on GitHub](https://github.com/pathogenwatch-oss/amr-search) |
        | Original Publication(s) | [PAARSNP: rapid genotypic resistance prediction for _Neisseria gonorrhoeae_](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7545138/) |
