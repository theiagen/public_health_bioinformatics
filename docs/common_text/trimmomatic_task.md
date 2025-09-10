??? task "`Trimmomatic`: Read Trimming (default)"
    Read proccessing is available via `Trimmomatic` by default.

    Trimmomatic trims low-quality regions of Illumina paired-end or single-end reads with a sliding window (with a default window size of 4, specified with `trim_window_size`), cutting once the average quality within the window falls below the `trim_quality_trim_score` (default of 20 for paired-end, 30 for single-end). The read is discarded if it is trimmed below `trim_minlen` (default of 75 for paired-end, 25 for single-end).

    !!! techdetails "`Trimmomatic` Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_trimmomatic.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_trimmomatic.wdl) |
        | Software Source Code | [Trimmomatic on GitHub](https://github.com/usadellab/Trimmomatic) |
        | Software Documentation | [Trimmomatic Website](http://www.usadellab.org/cms/?page=trimmomatic) |
        | Original Publication(s) | [Trimmomatic: a flexible trimmer for Illumina sequence data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4103590/) |
