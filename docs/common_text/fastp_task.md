<!-- if: read_qc_trim -->
??? task "`fastp`: Read Trimming (alternative)"
    To activate this task, set `read_processing` to `"fastp"`.

    `fastp` trims low-quality regions of Illumina paired-end or single-end reads with a sliding window (with a default window size of 4 (10 for TheiaEuk), specified with `trim_window_size`), cutting once the average quality within the window falls below the `trim_quality_trim_score` (default of 20 for paired-end, 30 for single-end). The read is discarded if it is trimmed below `trim_minlen` (default of 75 bases for paired-end, 25 for single-end).

    `fastp` also has additional default parameters and features that are not a part of `trimmomatic`'s default configuration. Please note `--disable_adapter_trimming` is explicitly needed in `fastp_args` to disable adapter trimming via `fastp`.

    ??? toggle "`fastp` default read-trimming parameters"
        | Parameter | Explanation |
        | --- | --- |
        | -g | enables polyG tail trimming |
        | -5 20 | enables read end-trimming |
        | -3 20 | enables read end-trimming |
        | --detect_adapter_for_pe | More sensitively detects adapters for trimming **only for paired-end reads** |
<!-- endif -->


<!-- if: theiaviral | metabuli -->
??? task "`fastp`: Read Trimming"
    `fastp` trims low-quality regions with a sliding window (with a default window size of 4, specified with `trim_window_size`), cutting once the average quality within the window falls below the `trim_quality_trim_score` (default of 20 for paired-end, 30 for single-end). The read is discarded if it is trimmed below `trim_minlen` (default of 15 bases).

    Adapter trimming is enabled by default and can be disabled by setting `fastp_trim_adapters` to "false".
<!-- endif -->

    Additional arguments can be passed using the `fastp_args` optional parameter. Please reference the [Fastp GitHub](https://github.com/OpenGene/fastp) for a comprehensive list of arguments.

    !!! techdetails "Trimmomatic and fastp Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_fastp.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_fastp.wdl) |
        | Software Source Code | [fastp on GitHub](https://github.com/OpenGene/fastp) |
        | Software Documentation | [fastp on GitHub](https://github.com/OpenGene/fastp) |
        | Original Publication(s) | [fastp: an ultra-fast all-in-one FASTQ preprocessor](https://academic.oup.com/bioinformatics/article/34/17/i884/5093234?login=false) |
