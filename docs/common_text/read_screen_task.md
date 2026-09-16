---
title: Task Fragment `read_screen`
fragment: true
---
??? task "`screen`: Total Raw Read Quantification and Genome Size Estimation"
    The screen task ensures the quantity of sequence data is sufficient to undertake genomic analysis. It uses [fastq-scan](https://github.com/rpetit3/fastq-scan) and bash commands for quantification of reads and base pairs, and [Mash](https://mash.readthedocs.io/en/latest/index.html) sketching to estimate the genome size and its coverage. At each step, the results are assessed relative to pass/fail criteria and thresholds that may be defined by optional user inputs. Samples are run through all threshold checks, regardless of failures, and the workflow will terminate after the screen task if any thresholds are not met:

    1. Total number of reads: A sample will fail the read screening task if its total number of reads is less than or equal to `min_reads`.
    2. The proportion of basepairs in the forward and reverse read files: A sample will fail the read screening if fewer than `min_proportion` basepairs are in either the reads1 or read2 files.
    3. Number of basepairs: A sample will fail the read screening if there are fewer than `min_basepairs` basepairs
    4. Estimated genome size:  A sample will fail the read screening if the estimated genome size is smaller than `min_genome_length` or bigger than `max_genome_length`.
    5. Estimated genome coverage: A sample will fail the read screening if the estimated genome coverage is less than the `min_coverage`.

<!-- if: theiacov -->
    Read screening is undertaken on both the raw and cleaned reads. Each screen has its own gate: set `skip_screen_raw` to `true` to skip the screen on the raw reads, and `skip_screen_clean` to `true` to skip the screen on the cleaned reads. The two are independent, so you can screen the raw reads only, the cleaned reads only, both, or neither.

<!-- endif -->

<!-- if: theiaprok|theiaeuk|theiaeukont -->
    Read screening is undertaken on both the raw and cleaned reads. The task may be skipped by setting the `skip_screen` variable to `true`.

    Default values vary between the PE, SE, and ONT workflows. The rationale for these default values can be found below. If two default values are shown, the first is for Illumina workflows and the second is for ONT.
<!-- endif -->

<!-- if: theiaviral -->
    Read screening is performed only on the cleaned reads. The task may be skipped by setting the `skip_screen` variable to `true`. Default values vary between the ONT and PE workflow. The rationale for these default values can be found below:

    ??? toggle "Default Thresholds and Rationales"
        | Variable  | Description | Default Value | Rationale |
        | --- | --- | --- | --- |
        | `estimated_genome_length` | | 12500 | This is an approximate of the median RNA virus length |
        | `min_reads` | A sample will fail the read screening task if its total number of reads is less than or equal to `min_reads` | 50 | Minimum number of base pairs for 10x coverage of the Hepatitis delta (of the _Deltavirus_ genus) virus divided by 300 (longest Illumina read length) |
        | `min_basepairs` | A sample will fail the read screening if there are fewer than `min_basepairs` basepairs | 15000 | Greater than 10x coverage of the Hepatitis delta (of the _Deltavirus_ genus) virus |
        | `min_genome_length` | A sample will fail the read screening if the estimated genome size is smaller than `min_genome_length` | 1500 |  Based on the Hepatitis delta (of the _Deltavirus_ genus) genome- the smallest viral genome as of 2024-04-11 (1,700 bp) |
        | `max_genome_length` | A sample will fail the read screening if the estimated genome size is larger than `max_genome_length` | 2673870 | Based on the _Pandoravirus salinus_ genome, the biggest viral genome, (2,673,870 bp) with 2 Mbp added |
        | `min_coverage` | A sample will fail the read screening if the estimated genome coverage is less than the `min_coverage` | 10 | A bare-minimum coverage for genome characterization. Higher coverage would be required for high-quality phylogenetics. |
        | `min_proportion` | A sample will fail the read screening if fewer than `min_proportion` basepairs are in either the reads1 or read2 files | 40 | Greater than 50% reads are in the read1 file; others are in the read2 file. (PE workflow only) |
<!-- endif -->

<!-- if: theiacov -->
    !!! warning "Breaking change: `skip_screen` has been replaced"
        The single `skip_screen` input has been split into `skip_screen_raw` and `skip_screen_clean`, and `skip_screen` is no longer a recognized input for TheiaCoV_Illumina_PE, TheiaCoV_Illumina_SE, or TheiaCoV_ONT.

        Existing workspace configurations and input JSONs that set `skip_screen` will **not** carry over. The value is silently dropped, both screens fall back to their `false` default, and samples that previously bypassed screening will now be screened and may terminate at the screen task. Replace `skip_screen` with `skip_screen_raw` and/or `skip_screen_clean` before rerunning. Other workflow series (TheiaProk, TheiaEuk, TheiaViral) are unaffected and keep the single `skip_screen` input.

    | Variable  | Default Value | Rationale |
    | --- | --- | --- |
    | `skip_screen_raw` | false | Set to true to skip the read screen on the raw reads. If set to true, `est_genome_length` is not estimated, so provide `genome_length` explicitly when `call_rasusa` is enabled |
    | `skip_screen_clean` | false | Set to true to skip the read screen on the cleaned reads |
    | `min_reads` | 57 | Calculated from the minimum number of base pairs required for 10x coverage of the Hepatitis delta (of the _Deltavirus_ genus) genome, the smallest known viral genome as of 2024-04-11 (1,700 bp), divided by 300 (the longest Illumina read length) |
    | `min_basepairs` | 17000 | Should be greater than 10x coverage of Hepatitis delta (of the _Deltavirus_ genus), the smallest known viral genome (1,700 bp) |
    | `min_genome_length` | 1700 | Based on the Hepatitis delta (of the _Deltavirus_ genus) genome, the smallest viral genome as of 2024-04-11 (1,700 bp) |
    | `max_genome_length` | 2673870 | Based on the _Pandoravirus salinus_ genome, the largest known viral genome (2,673,870 bp), plus an additional 200 kbp to cater for potential extra genomic material |
    | `min_coverage` | 10 | A bare-minimum average per base coverage across the genome required for genome characterization. Higher coverage would be required for high-quality phylogenetics. |
    | `min_proportion` | 40 | Neither read1 nor read2 files should have less than 40% of the total number of reads. For paired-end data only. |
<!-- endif -->

<!-- if: theiaprok -->
    | Variable  | Default Value | Rationale |
    | --- | --- | --- |
    | `skip_screen` | false | Set to true to skip the read screen from running |
    | `min_reads` | 7472 or 5000 | Calculated from the minimum number of base pairs required for 20x coverage of the Nasuia deltocephalinicola genome, the smallest known bacterial genome as of 2019-08-07 (112,091 bp), divided by 300 (the longest Illumina read length) or 5000 (estimate of ONT read length) |
    | `min_basepairs` | 2241820 | Should be greater than 20x coverage of Nasuia deltocephalinicola, the smallest known bacterial genome (112,091 bp) |
    | `min_genome_length` | 100000 | Based on the Nasuia deltocephalinicola genome, the smallest known bacterial genome (112,091 bp) |
    | `max_genome_length` | 18040666 | Based on the Minicystis rosea genome, the largest known bacterial genome (16,040,666 bp), plus an additional 2 Mbp to cater for potential extra genomic material |
    | `min_coverage` | 10 or 5 | A bare-minimum average per base coverage across the genome required for genome characterization. Higher coverage would be required for high-quality phylogenetics. |
    | `min_proportion` | 40 | Neither read1 nor read2 files should have less than 40% of the total number of reads. For paired-end data only. |
<!-- endif -->

<!-- if: theiaeuk -->
    | Variable  | Default Value | Rationale |
    | --- | --- | --- |
    | `skip_screen` | false | Set to true to skip the read screen from running. If you set this value to true, please provide a value for the theiaeuk_illumina_pe `genome_length` optional input, OR set the theiaeuk_illumina_pe `call_rasusa` optional input to false. Otherwise Rasusa will attempt to downsample to an expected genome size of 0 bp, and the workflow will fail. |
    | `min_reads` | 30000 | Calculated from the minimum number of base pairs required for 20x coverage of the _Hansenula polymorpha_ genome, the smallest fungal genome as of 2015-04-02 (8.97 Mbp), divided by 300 (the longest Illumina read length) |
    | `min_basepairs` | 45000000 | Should be greater than 10x coverage of _Hansenula polymorpha_, the smallest fungal genome as of 2015-04-02 (8.97 Mbp)  |
    | `min_genome_length` | 9000000 | Based on the _Hansenula polymorpha_  genome - the smallest fungal genome as of 2015-04-02 (8.97 Mbp) |
    | `max_genome_length` | 178000000 | Based on the _Cenococcum geophilum_  genome, the largest pathogenic fungal genome (177.57 Mbp), plus an additional 2 Mbp to cater for potential extra genomic material |
    | `min_coverage` | 10 | A bare-minimum average per base coverage across the genome required for genome characterization. Higher coverage would be required for high-quality phylogenetics.|
    | `min_proportion` | 40 | Neither read1 nor read2 files should have less than 40% of the total number of reads. For paired-end data only. |
<!-- endif -->

<!-- if: theiaeukont -->
    | Variable  | Default Value | Rationale |
    | --- | --- | --- |
    | `skip_screen` | false | Set to true to skip the read screen from running. |
    | `genome_length` | 50000000 | An average genome length to avoid erroneous estimated genome length calculations. It is recommended to provide the expected genome length if available. |
    | `min_reads` | 5000 | Calculated from the minimum number of base pairs required for 20x coverage of the _Hansenula polymorpha_ genome, the smallest fungal genome as of 2015-04-02 (8.97 Mbp), divided by 1800 |
    | `min_basepairs` | 45000000 | Should be greater than 10x coverage of _Hansenula polymorpha_, the smallest fungal genome as of 2015-04-02 (8.97 Mbp)  |
    | `min_genome_length` | 9000000 | Based on the _Hansenula polymorpha_  genome - the smallest fungal genome as of 2015-04-02 (8.97 Mbp) |
    | `max_genome_length` | 178000000 | Based on the _Cenococcum geophilum_  genome, the largest pathogenic fungal genome (177.57 Mbp), plus an additional 2 Mbp to cater for potential extra genomic material |
    | `min_coverage` | 5 | A bare-minimum average per base coverage across the genome required for genome characterization. Higher coverage would be required for high-quality phylogenetics.|
<!-- endif -->

    !!! techdetails "Screen Technical Details"

<!-- if: theiacov|theiaprok|theiaeuk|theiaeukont -->
        There is a single task for read screening. The screen task is run twice, once for raw reads and once for clean reads.
<!-- endif -->

        |  | Links |
        | --- | --- |
        | Task | [task_screen.wdl (PE sub-task)](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/comparisons/task_screen.wdl#L3)<br>[task_screen.wdl (SE sub-task)](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/comparisons/task_screen.wdl#L147) |
