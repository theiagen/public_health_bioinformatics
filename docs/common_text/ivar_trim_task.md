??? task "`ivar_trim`: Primer Trimming (optional)"
    To deactivate this task, set `trim_primers` to `false`.

    Using the user-provided (or, more rarely, a [_organism-specific parameters_-determined](./theiacov.md#org-specific)) `primer_bed` file, iVar soft-clips primer sequences from an aligned and sorted BAM file and then trims the reads based on a quality threshold of 20 using a sliding window approach. If the resulting read is greater than 30 bp, the read is written to a a new BAM file consisting of only trimmed reads (or reads that did not have a primer identified).

    !!! techdetails "iVar Trim Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_ivar_primer_trim.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_ivar_primer_trim.wdl) |
        | Software Source Code | [iVar on GitHub](https://andersen-lab.github.io/ivar/html/) |
        | Software Documentation | [iVar on GitHub](https://andersen-lab.github.io/ivar/html/) |
        | Original Publication(s) | [An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar](http://dx.doi.org/10.1186/s13059-018-1618-7) |
