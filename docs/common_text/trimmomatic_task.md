---
title: Task Fragment `trimmomatic`
fragment: true
---
??? task "`Trimmomatic`: Read Trimming (default)"
    Read processing is available via Trimmomatic by default.

    Trimmomatic trims low-quality regions of Illumina paired-end or single-end reads with a sliding window (with a default window size of 4, specified with `trim_window_size`), cutting once the average quality within the window falls below the `trimmomatic_window_quality` (default of 30 for both paired-end and single-end). The read is discarded if it is trimmed below `trimmomatic_min_length` (default of 75 for paired-end, 25 for single-end).

    Adapter trimming with Trimmomatic is disabled by default. It can be enabled by setting `trimmomatic_trim_adapters=true`. When enabled, Trimmomatic uses its default adapter sequences, [`TruSeq3-PE-2.fa`](https://github.com/usadellab/Trimmomatic/tree/main/adapters) (for paired-end) or [`TruSeq3-SE.fa`](https://github.com/usadellab/Trimmomatic/tree/main/adapters) (for single-end) with default adapter clipping settings equivalent to:

    ```
    ILLUMINACLIP:<adapter_fasta>:2:30:10

    2 = seed mismatches
    30 = palindrome clip threshold
    10 = simple clip threshold
    ```

    Users can optionally provide a custom adapter file or modify adapter trimming parameters using the `trimmomatic_adapter_fasta` and `trimmomatic_adapter_trim_args` respectively. See the [Trimmomatic adapter documentation](https://github.com/usadellab/Trimmomatic?tab=readme-ov-file#the-adapter-fasta-files) for more details. The `trimmomatic_adapter_fasta` parameter should just include the path to your fasta file. The `trimmomatic_adapter_trim_args` parameter should only contain the colon-delimited values that come after the adapter fasta file in the `ILLUMINACLIP` argument. Example usage:

    ```
    trimmomatic_adapter_fasta="path/to/my_custom_adapters.fa"
    trimmomatic_adapter_trim_args="2:30:10"
    ```

    Trimmomatic's automatic detection of the phred quality score encoding is not always reliable, particularly for older or SRA-Lite-formatted FASTQ files. PHB explicitly sets `"-phred33"` as the default via `trimmomatic_quality_encoding`, since the vast majority of modern FASTQ files use `"-phred33"` encoding. The accepted values for `trimmomatic_quality_encoding` are `"-phred33"`, `"-phred64"`, and `""` (an empty string, which restores Trimmomatic's auto-detection). Note that the leading dash is required and `"-phred64"` should only be used for genuinely pre-2011 Illumina data.

    ```
    trimmomatic_quality_encoding="-phred33"
    ```

    For more advanced configurations, there are options to override the default trimming parameters via `trimmomatic_override_args`. Note that when using `trimmomatic_override_args`, the user is responsible for specifying all desired trimming steps and their order, as the default trimming steps will be ignored. Adapter trimming and the quality score encoding are configured separately and are not affected. See the [Trimmomatic documentation](https://github.com/usadellab/Trimmomatic?tab=readme-ov-file#step-options) for more details on available trimming steps and their parameters.

    Advanced Configuration Example Usage:
    ```
    trimmomatic_override_args="LEADING:30 TRAILING:30 AVGQUAL:36 SLIDINGWINDOW:4:33 MINLEN:75"
    ```

    !!! techdetails "Trimmomatic Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_trimmomatic.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_trimmomatic.wdl) |
        | Software Source Code | [Trimmomatic on GitHub](https://github.com/usadellab/Trimmomatic) |
        | Software Documentation | [Trimmomatic Website](http://www.usadellab.org/cms/?page=trimmomatic) |
        | Original Publication(s) | [Trimmomatic: a flexible trimmer for Illumina sequence data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4103590/) |
