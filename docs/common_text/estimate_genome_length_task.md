??? task "`datasets_genome_length`"
    The `datasets_genome_length` task uses [`NCBI Datasets`](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/taxonomy/) to acquire genome length metadata for an inputted taxon and retrieve a top reference accession. This task generates a summary file of all successful hits to the input `taxon`, which includes each genome's accession number, completeness status, genome length, source, and other relevant metadata. The task will then calculate the average expected genome length in basepairs for the input `taxon`.

    ??? dna "`taxon` input parameter"
        This parameter accepts either a NCBI taxon ID (e.g. `11292`) or an organism name (e.g. `Lyssavirus rabies`).

    !!! techdetails "NCBI Datasets Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_identify_taxon_id.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/task_identify_taxon_id.wdl) |
        | Software Source Code | [NCBI Datasets on GitHub](https://github.com/ncbi/datasets) |
        | Software Documentation | [NCBI Datasets Documentation on NCBI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/) |
        | Original Publication(s) | [Exploring and retrieving sequence and metadata for species across the tree of life with NCBI Datasets](https://doi.org/10.1038/s41597-024-03571-y) |
