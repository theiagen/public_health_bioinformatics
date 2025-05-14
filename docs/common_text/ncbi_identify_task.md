<!-- if: theiaviral -->
??? task "`ncbi_identify`"

    The `ncbi_identify` task uses [`NCBI Datasets`](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/taxonomy/) to search the NCBI Viral Genome Database and acquire taxonomic metadata from a user's inputted taxonomy and desired taxonomic rank. This task will always return a taxon ID, name, and rank, and it facilitates multiple downstream functions, including **read classification and targeted read extraction**.

    ??? dna "`taxon`"
        This parameter accepts either a NCBI taxon ID (e.g. `11292`) or an organism name (e.g. `*Lyssavirus rabies*`).

    ??? dna "`rank` a.k.a `read_extraction_rank`"
        Valid options include: `"species"`, `"genus"`, `"family"`, `"order"`, `"class"`, `"phylum"`, `"kingdom"`, or `"domain"`. By default it is set to `"family"`. This parameter filters metadata to report information only at the taxonomic `rank` specified by the user, regardless of the taxonomic level implied by the original input `taxon`.

    ???+ warning "Important"
        - The `rank` parameter must specify a taxonomic level that is equal to or broader than the input organism's taxonomic level.

        **Examples:**

        - If your input `taxon` is `"*Lyssavirus rabies*"` (species level) with `rank` set to `"family"`, the task will return the family-level information: taxon ID for Rhabdoviridae (11270), name "Rhabdoviridae", and rank "family".
        - If your input `taxon` is `"*Lyssavirus*"` (genus level) with `rank` set to `"species"`, the task will fail because it cannot determine species-level information from a genus-level input.
<!-- endif -->

    !!! techdetails "NCBI Datasets Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_identify_taxon_id.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/task_identify_taxon_id.wdl) |
        | Software Source Code | [NCBI Datasets on GitHub](https://github.com/ncbi/datasets) |
        | Software Documentation | [NCBI Datasets Documentation on NCBI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/) |
        | Original Publication(s) | [Exploring and retrieving sequence and metadata for species across the tree of life with NCBI Datasets](https://doi.org/10.1038/s41597-024-03571-y) |
