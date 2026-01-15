??? task "`Vibecheck`: O1 _Vibrio  cholerae_ Lineage Classification ==_for Illumina PE only_=="
    The `Vibecheck` task classifies O1 _V. cholerae_ sequences into canonical lineages (T1-T17) using variant frequency demixing. The O1 designation is determined through the use of the SRST2 task.

    Vibecheck works by aligning the reads to an O1 _V. cholerae_ reference genome, calling variants from the alignment, and estimating lineage abundances using [Freyja](./freyja.md) by using a database built from canonical SNPs that define the known lineages.

    !!! techdetails "Vibecheck Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_vibecheck_vibrio.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/vibrio/task_vibecheck_vibrio.wdl) |
        | Software Source Code | [Vibecheck on GitHub](https://github.com/CholGen/Vibecheck) |
        | Software Documentation | [Vibecheck on GitHub](https://github.com/CholGen/Vibecheck) |
