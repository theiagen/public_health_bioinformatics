---
title: Task Fragment `ivar_trim`
fragment: true
---
<!-- if: freyja -->
??? task "`iVar trim`: Primer Trimming"
    The optional input, `keep_noprimer_reads`, does not have to be modified.
<!-- endif -->
<!-- if: theiacov -->
??? task "`iVar trim`: Primer Trimming (optional)"
    To deactivate this task, set `trim_primers` to `false`.
<!-- endif -->
<!-- if: theiaviral -->
??? task "`iVar trim`: Primer Trimming (optional)"
    To activate this task, provide a `primer_bed` file containing (0-based index) primer coordinates in BED format.
<!-- endif -->

    Using the user-provided (or a [_organism-specific parameters_-determined](./theiacov.md#org-specific)) `primer_bed` file, iVar soft-clips primer sequences from an aligned and sorted BAM file. 
    
    iVar will trim any reads that start or end within the (0-based index) coordinates provided in the BED file. It does not take the sequence of bases itself into account. This allows iVar to accurately trim primer sequences despite potential mismatches between sequencing reads and primer sequences in the aligned region.

    Following the trimming of primer sequences, iVar then trims the reads based on a quality threshold of 20 using a sliding window approach (default: 4). If the average base quality drops below the threshold, the remainder of the read is soft-clipped. Reads exceeding the minimum length (default: 30) after trimming are retained in the output BAM file.

    !!! techdetails "iVar Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_ivar_primer_trim.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_ivar_primer_trim.wdl) |
        | Software Source Code | [iVar on GitHub](https://andersen-lab.github.io/ivar/html/) |
        | Software Documentation | [iVar on GitHub](https://andersen-lab.github.io/ivar/html/) |
        | Original Publication(s) | [An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar](http://dx.doi.org/10.1186/s13059-018-1618-7) |
