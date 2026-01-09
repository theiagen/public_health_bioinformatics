??? task "`ivar_trim`: Primer Trimming (optional)"

<!-- if: theiacov -->
    To deactivate this task, set `trim_primers` to `false`.

    Using the user-provided (or, more rarely, a [_organism-specific parameters_-determined](./theiacov.md#org-specific)) `primer_bed` file, iVar soft-clips primer sequences from an aligned and sorted BAM file and then trims the reads based on a quality threshold of 20 using a sliding window approach. If the resulting read is greater than 30 bp, the read is written to a a new BAM file consisting of only trimmed reads (or reads that did not have a primer identified).
<!-- endif -->

<!-- if: theiaviral -->
    To activate this task, provide a `primer_bed` file containing (0-based) primer coordinates in BED format.

    First, iVar soft-clips primer sequences from the aligned reads in a sorted BAM file using the coordinates provided by a user-supplied `primer_bed` file. iVar will trim any reads that start or end within the (0-based) coordinates provided in the BED file. It does not take the sequence of bases itself into account. This allows iVar to accurately trim primer sequences despite potential mismatches between sequencing reads and primer sequences in the aligned region.

    Following the trimming of primer sequences, iVar then trims the reads based on a quality threshold of 20 using a sliding window approach (default: 4). If the average base quality drops below the threshold, the remainder of the read is soft-clipped. Reads exceeding the minimum length (default: 30) after trimming are retained in the output BAM file.
<!-- endif -->

    !!! techdetails "iVar Trim Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_ivar_primer_trim.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_ivar_primer_trim.wdl) |
        | Software Source Code | [iVar on GitHub](https://andersen-lab.github.io/ivar/html/) |
        | Software Documentation | [iVar on GitHub](https://andersen-lab.github.io/ivar/html/) |
        | Original Publication(s) | [An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar](http://dx.doi.org/10.1186/s13059-018-1618-7) |
