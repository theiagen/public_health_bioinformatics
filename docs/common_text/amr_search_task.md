<!-- if: amr_search -->
??? task "`amr_search`: Antimicrobial resistance profiling"
<!-- endif -->

<!-- if: theiaprok -->
??? task "`amr_search`: _Neisseria gonorrhoeae_ antimicrobial resistance profiling"

    This task performs *in silico* antimicrobial resistance (AMR) profiling for *Neisseria gonorrhoeae* using **AMRsearch**, the primary tool used by [Pathogenwatch](https://pathogen.watch/) to genotype and infer antimicrobial resistance (AMR) phenotypes from assembled microbial genomes.

    The AMR search is conducted when *Neisseria gonorrhoeae* is identified as the taxon in *TheiaProk* workflows. The default database for *N. gonorrhoeae* is **485**.

<!-- endif -->

<!-- if: theiaeuk -->
??? task "`amr_search`: Antimicrobial resistance profiling (Optional)"

    Set the `run_amr_search` parameter to `true` to enable this task.

    The AMR Search module will report only specific SNP changes in the AMR genes, like FUR1 F211L. For a complete list of mutations being queried, please refer to [`amr_search`'s documentation](https://github.com/pathogenwatch-oss/amr-search). 
<!-- endif -->

<!-- if: amr_search|theiaeuk -->
    This task performs *in silico* antimicrobial resistance (AMR) profiling for supported species using **AMRsearch**, the primary tool used by [Pathogenwatch](https://pathogen.watch/) to genotype and infer antimicrobial resistance (AMR) phenotypes from assembled microbial genomes.
<!-- endif -->

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