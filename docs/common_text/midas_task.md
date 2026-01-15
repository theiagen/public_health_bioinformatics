??? task "`MIDAS`: Read Identification (optional)"
    To activate this task, set `call_midas` to `true`.

    The `MIDAS` task is for the identification of reads to detect contamination with non-target taxa.

    The MIDAS tool was originally designed for metagenomic sequencing data but has been co-opted for use with bacterial isolate WGS methods. It can be used to detect contamination present in raw sequencing data by estimating bacterial species abundance in bacterial isolate WGS data. If a secondary genus is detected above a relative frequency of 0.01 (1%), then the sample should fail QC and be investigated further for potential contamination.

    This task is similar to those used in commercial software, BioNumerics, for estimating secondary species abundance.

    ??? toggle "How are the MIDAS output columns determined?"
        
        Example MIDAS report in the `midas_report` column:
        
        | species_id | count_reads | coverage | relative_abundance |
        | --- | --- | --- | --- |
        | Salmonella_enterica_58156 | 3309 | 89.88006645 | 0.855888033 |
        | Salmonella_enterica_58266 | 501 | 11.60606061 | 0.110519371 |
        | Salmonella_enterica_53987 | 99 | 2.232896237 | 0.021262881 |
        | Citrobacter_youngae_61659 | 46 | 0.995216227 | 0.009477003 |
        | Escherichia_coli_58110 | 5 | 0.123668877 | 0.001177644 |
        
        MIDAS report column descriptions:
        
        - species_id: species identifier
        - count_reads: number of reads mapped to marker genes
        - coverage: estimated genome-coverage (i.e. read-depth) of species in metagenome
        - relative_abundance: estimated relative abundance of species in metagenome
        
        The value in the `midas_primary_genus` column is derived by ordering the rows in order of "relative_abundance" and identifying the genus of top species in the "species_id" column (Salmonella). The value in the `midas_secondary_genus` column is derived from the genus of the second-most prevalent genus in the "species_id" column (Citrobacter). The `midas_secondary_genus_abundance` column is the "relative_abundance" of the second-most prevalent genus (0.009477003). The `midas_secondary_genus_coverage` is the "coverage" of the second-most prevalent genus (0.995216227).

    **MIDAS Reference Database Overview**

    The **MIDAS reference database** is a comprehensive tool for identifying bacterial species in metagenomic and bacterial isolate WGS data. It includes several layers of genomic data, helping detect species abundance and potential contaminants.

    !!! dna "Key Components of the MIDAS Database"
        1. **Species Groups**: 
            - MIDAS clusters bacterial genomes based on 96.5% sequence identity, forming over 5,950 species groups from 31,007 genomes. These groups align with the gold-standard species definition (95% ANI), ensuring highly accurate species identification.

        2. **Genomic Data Structure**:
            - _Marker Genes_: Contains 15 universal single-copy genes used to estimate species abundance.
            - _Representative Genome_: Each species group has a selected representative genome, which minimizes genetic variation and aids in accurate SNP identification.
            - _Pan-genome_: The database includes clusters of non-redundant genes, with options for multi-level clustering (e.g., 99%, 95%, 90% identity), enabling MIDAS to identify gene content within strains at various clustering thresholds.

        3. **Taxonomic Annotation**: 
            - Genomes are annotated based on consensus Latin names. Discrepancies in name assignments may occur due to factors like unclassified genomes or genus-level ambiguities.

    ---

    **Using the Default MIDAS Database**

    TheiaProk and TheiaEuk use the pre-loaded MIDAS database in Terra (see input table for current version) by default for bacterial species detection in metagenomic data, requiring no additional setup.

    !!! tip "Create a Custom MIDAS Database"
        Users can also build their own custom MIDAS database if they want to include specific genomes or configurations. This custom database can replace the default MIDAS database used in Terra. To build a custom MIDAS database, follow the [MIDAS GitHub guide on building a custom database](https://github.com/snayfach/MIDAS/blob/master/docs/build_db.md). Once the database is built, users can upload it to a Google Cloud Storage bucket or Terra workkspace and provide the link to the database in the `midas_db` input variable.

    ---

    !!! techdetails "MIDAS Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_midas.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/contamination/task_midas.wdl) |
        | Software Source Code | [MIDAS on GitHub](https://github.com/snayfach/MIDAS)
        | Software Documentation | [MIDAS on GitHub](https://github.com/snayfach/MIDAS)
        | Original Publication(s) | [An integrated metagenomics pipeline for strain profiling reveals novel patterns of bacterial transmission and biogeography](https://pubmed.ncbi.nlm.nih.gov/27803195) |