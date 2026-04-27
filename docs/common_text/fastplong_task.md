??? task "`fastplong`: ONT Read Trimming"
    `fastplong` trims low-quality regions with a sliding window (with a default window size of 4, specified with `fastplong_window_size`), cutting once the average quality within the window falls below the `fastplong_quality_trim_score` (default of 20). The read is discarded if it is trimmed below `fastplong_min_length` (default of 15 bases). These trimming options are conducted according to a sliding window, but the directionality of this window can be specified by setting `cut_front` (5' to 3') or `cut_tail` (3' to 5') to "true".

    Adapter trimming is enabled by default and can be disabled by setting `fastplong_trim_adapters` to "false". Automatic adapter detection is enabled by default, though a FASTA, a string of start, or a string of end adapters can be specified with `fastplong_adapter_fasta`, `fastplong_start_adapter`, or `fastplong_end_adapter` inputs respectively.

    Additional arguments can be passed using the `fastplong_args` optional parameter. Please reference the [Fastp GitHub](https://github.com/OpenGene/fastplong) for a comprehensive list of arguments.

    !!! techdetails "fastplong Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_fastp.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_fastplong.wdl) |
        | Software Source Code | [fastp on GitHub](https://github.com/OpenGene/fastplong) |
        | Software Documentation | [fastp on GitHub](https://github.com/OpenGene/fastplong) |
        | Original Publication(s) | [Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp](https://doi.org/10.1002/imt2.107) |
