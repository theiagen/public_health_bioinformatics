??? toggle "`HRRT`: Human Host Sequence Removal"

    All reads of human origin **are removed**, including their mates, by using NCBI's [**human read removal tool (HRRT)**](https://github.com/ncbi/sra-human-scrubber). 

    HRRT is based on the [SRA Taxonomy Analysis Tool](https://doi.org/10.1186/s13059-021-02490-0) and employs a k-mer database constructed of k-mers from Eukaryota derived from all human RefSeq records with any k-mers found in non-Eukaryota RefSeq records subtracted from the database.

    !!! techdetails "NCBI-Scrub Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_ncbi_scrub.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_ncbi_scrub.wdl) |
        | Software Source Code | [NCBI Scrub on GitHub](https://github.com/ncbi/sra-human-scrubber) |
        | Software Documentation | <https://github.com/ncbi/sra-human-scrubber/blob/master/README.md> |
