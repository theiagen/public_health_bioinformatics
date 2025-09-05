
??? task "`Pilon`: Assembly Polishing"
    `Pilon` is a tool that uses read alignments to correct errors in an assembly.

<!-- if: digger -->
    The `bwa`-generated alignment of the read data to the assembly is used to identify inconsistences between the reads and the assembly in order to correct them. `Pilon` will attempt to fix individual base errors and small indels using the read data. This can improve the overall quality of the assembly, especially when the assembler has made mistakes due to sequencing errors or low coverage regions.

    The default parameters were set to mimic the parameters used by [Shovill](https://github.com/tseemann/shovill): `--fix bases --minq 60 --minqual 3 --mindepth 0.25`. These can be modified by the user.
<!-- endif -->
<!-- if: theiameta -->
    It is used to polish the assembly produced by metaSPAdes. The input to Pilon is the sorted BAM file produced by `samtools`, and the original draft assembly produced by `metaspades`.
<!-- endif -->

    !!! techdetails "Pilon Technical Details"
        | | Links |
        |---|---|
        | Task | [task_pilon.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_pilon.wdl) |
        | Software Source Code | [Pilon on GitHub](https://github.com/broadinstitute/pilon) |
        | Software Documentation | [Pilon Wiki](https://github.com/broadinstitute/pilon/wiki) |
        | Original Publication(s) | [Pilon: An Integrated Tool for Comprehensive Microbial Variant Detection and Genome Assembly Improvement](https://doi.org/10.1371/journal.pone.0112963) |
