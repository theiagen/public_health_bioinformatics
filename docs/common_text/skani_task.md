??? task "`skani`"

    The `skani` task is used to identify and select the most closely related reference genome to the *de novo* assembly. Skani uses an approximate mapping method without base-level alignment to calculate average nucleotide identity (ANI). It is magnitudes faster than BLAST-based methods and almost as accurate.

<!-- if: theiaviral -->
    By default, the reference genome is selected from a database of approximately 200,000 complete viral genomes. This database was constructed with the following methodology:
    
    1. Extracting all [complete NCBI viral genomes](https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/AllNuclMetadata/), excluding RefSeq accessions (redundancy), SARS-CoV-2 accessions, and segmented families (Orthomyxoviridae, Hantaviridae, Arenaviridae, and Phenuiviridae)
    
    2. Adding complete RefSeq segmented viral assembly accessions, which represent segments as individual contigs within the FASTA

    3. Adding one SARS-CoV-2 genome for each [major pangolin lineages](https://github.com/cov-lineages/lineages-website/raw/refs/heads/master/_data/lineage_data.full.json)
    
<!-- endif -->
    !!! techdetails "Skani Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_skani.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/task_skani.wdl) |
        | Software Source Code | [Skani on GitHub](https://github.com/bluenote-1577/skani) |
        | Software Documentation | [Skani Documentation](https://github.com/bluenote-1577/skani/blob/main/README.md) |
        | Original Publication(s) | [Skani Paper](https://doi.org/10.1038/s41592-023-02018-3) |