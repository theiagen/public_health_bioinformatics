??? task "NCBI Datasets"
    ##### NCBI Datasets

    The [`NCBI Datasets`](https://www.ncbi.nlm.nih.gov/datasets/) task downloads specified assemblies from NCBI using either the [virus](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/virus-genome/) or [genome](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/genome/) (for all other genome types) package as appropriate.

<!-- if: snippy -->
    !!! dna "`include_gbff` behavior"
        If `include_gbff` is set to `true`, the gbff file will be used as the reference for `Snippy_Variants` and `Snippy_Tree`. If `include_gbff` is set to `false`, the fasta file will be used as the reference for `Snippy_Variants` and `Snippy_Tree`. Tree topology should not differ, though annotations may.
<!-- endif -->

    !!! techdetails "NCBI Datasets Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_ncbi_datasets.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/data_import/task_ncbi_datasets.wdl) |
        | Software Source Code | [NCBI Datasets on GitHub](https://github.com/ncbi/datasets) |
        | Software Documentation | [NCBI Datasets Documentation on NCBI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/) |
        | Original Publication(s) | [Exploring and retrieving sequence and metadata for species across the tree of life with NCBI Datasets](https://doi.org/10.1038/s41597-024-03571-y) |