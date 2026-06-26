---
title: Task Fragment `ete4_identify`
fragment: true
---
??? task "`ete4_identify`: Taxonomic Identification"
    The ETE Toolkit uses ete4 to parse the NCBI taxonomy hierarchy from a user's inputted taxonomy and desired taxonomic rank. This task returns a taxon ID, name, and rank, which facilitates downstream functions, including **read classification, targeted read extraction, and genomic characterization modules**.

    ??? dna "`taxon` input parameter"
        This parameter accepts either an NCBI taxon ID (e.g. `11292`) or an organism name (e.g. `Lyssavirus rabies`).

    ??? dna "`rank` a.k.a `read_extraction_rank` input parameter"
        Valid options include: `"species"`, `"genus"`, `"family"`, `"order"`, `"class"`, `"phylum"`, `"kingdom"`, or `"domain"`. By default it is set to `"family"`. This parameter filters metadata to report information only at the taxonomic `rank` specified by the user, regardless of the taxonomic rank implied by the original input `taxon`.

    ???+ warning "Important"
        - The `rank` parameter **must** specify a taxonomic rank that is ==equal to or above== the input taxon's taxonomic rank.

        **Examples:**

        - If your input `taxon` is `Lyssavirus rabies` (species level) with `rank` set to `family`, the task will return information for the family of `Lyssavirus rabies`: taxon ID for Rhabdoviridae (11270), name "Rhabdoviridae", and rank "family".
        - If your input `taxon` is `Lyssavirus` (genus level) with `rank` set to `species`, the task will fail because it cannot determine species information from an inputted genus.

    !!! techdetails "ETE Toolkit Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_ete4_taxon_id.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/task_ete4_taxon_id.wdl) |
        | Software Source Code | [ETE Toolkit on GitHub](https://github.com/etetoolkit/ete) |
        | Software Documentation | [ETE Toolkit Documentation](https://etetoolkit.github.io/ete/) |
        | Original Publication(s) | [ETE 3: Reconstruction, analysis and visualization of phylogenomic data](https://doi.org/10.1093/molbev/msw046) |
