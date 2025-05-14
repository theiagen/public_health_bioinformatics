??? task "`skani`"

<!-- if: theiaviral -->
    The `skani` task is used to identify and select the most closely related reference genome to the input assembly generated from the `spades` or `megahit` tasks. This reference genome is selected from a comprehensive database of over 270,000 complete viral genomes. Skani uses an approximate mapping method without base-level alignment to get ANI. It is magnitudes faster than BLAST-based methods and almost as accurate.
<!-- endif -->
    !!! techdetails "Skani Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_skani.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/task_skani.wdl) |
        | Software Source Code | [Skani on GitHub](https://github.com/bluenote-1577/skani) |
        | Software Documentation | [Skani Documentation](https://github.com/bluenote-1577/skani/blob/main/README.md) |
        | Original Publication(s) | [Skani Paper](https://doi.org/10.1038/s41592-023-02018-3) |