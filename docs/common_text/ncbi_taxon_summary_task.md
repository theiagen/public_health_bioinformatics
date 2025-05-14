<!-- if: theiaviral -->
??? task "`ncbi_taxon_summary`"

    The `ncbi_taxon_summary` task uses [`NCBI Datasets`](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/virus-genome/) to search the NCBI Viral Genome Database and acquire metadata based on a user's taxonomic input. This task generates a comprehensive summary file of all successful hits to the input `taxon`, which includes each taxon's accession number, completeness status, genome length, source, and other relevant metadata. Based on this summary, the task also calculates the average expected genome size for the input `taxon`.

    ??? dna "`taxon`"
        This parameter accepts either a NCBI taxon ID (e.g. `11292`) or an organism name (e.g. `Lyssavirus rabies`).

    ???+ warning "Important"
        - The implied taxonomic level of the input `taxon` directly affects the number of taxa reported in the summary file and the average genome size calculation.

        **Examples:**

        - If your input `taxon` is `"*Lyssavirus rabies*"` or taxon ID `11292` (species level), the task will retrieve genomes only for that specific species, and return the genome length for that species.
        - If your input `taxon` is `"Rhabdoviridae"` or taxon ID `11270` (family level), the task will retrieve genomes for all species within that family and calculate an average genome length across all retrieved species.

        This flexibility allows users to run the workflow without requiring precise knowledge of their organism's genome length. The average genome length calculation is only used to estimate coverage levels for downsampling and minor read QC steps. The average genome length is used only to estimate coverage levels for downsampling and to guide minor read quality control steps. It does not significantly affect the quality of the consensus genome or the choice of reference sequences.
<!-- endif -->

    !!! techdetails "NCBI Datasets Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_ncbi_datasets.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/data_import/task_ncbi_datasets.wdl) |
        | Software Source Code | [NCBI Datasets on GitHub](https://github.com/ncbi/datasets) |
        | Software Documentation | [NCBI Datasets Documentation on NCBI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/) |
        | Original Publication(s) | [Exploring and retrieving sequence and metadata for species across the tree of life with NCBI Datasets](https://doi.org/10.1038/s41597-024-03571-y) |