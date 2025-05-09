
??? task "`screen`: Total Raw Read Quantification and Genome Size Estimation"

    The [`screen`](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/comparisons/task_screen.wdl) task ensures the quantity of sequence data is sufficient to undertake genomic analysis. It uses [`fastq-scan`](https://github.com/rpetit3/fastq-scan) and bash commands for quantification of reads and base pairs, and [mash](https://mash.readthedocs.io/en/latest/index.html) sketching to estimate the genome size and its coverage. At each step, the results are assessed relative to pass/fail criteria and thresholds that may be defined by optional user inputs. Samples are run through all threshold checks, regardless of failures, and the workflow will terminate after the `screen` task if any thresholds are not met:

    1. Total number of reads: A sample will fail the read screening task if its total number of reads is less than or equal to `min_reads`.
    2. The proportion of basepairs reads in the forward and reverse read files: A sample will fail the read screening if fewer than `min_proportion` basepairs are in either the reads1 or read2 files.
    3. Number of basepairs: A sample will fail the read screening if there are fewer than `min_basepairs` basepairs
    4. Estimated genome size:  A sample will fail the read screening if the estimated genome size is smaller than `min_genome_size` or bigger than `max_genome_size`.
    5. Estimated genome coverage: A sample will fail the read screening if the estimated genome coverage is less than the `min_coverage`.

    Read screening is undertaken on both the raw and cleaned reads. The task may be skipped by setting the `skip_screen` variable to true.

    Default values vary between the PE, SE, and ONT workflows. The rationale for these default values can be found below. If two default values are shown, the first is for Illumina workflows and the second is for ONT.

<!-- if: theiacov -->
    | Variable  | Default Value | Rationale |
    | --- | --- | --- |
    | `skip_screen` | false | Set to true to skip the read screen from running |
    | `min_reads` | 57 | Calculated from the minimum number of base pairs required for 10x coverage of the Hepatitis delta (of the _Deltavirus_ genus) genome, the samllest known viral genome as of 2024-04-11 (1,700 bp), divided by 300 (the longest Illumina read length) |
    | `min_basepairs` | 17000 | Should be greater than 10x coverage of Hepatitis delta (of the _Deltavirus_ genus), the smallest known viral genome (1,700 bp) |
    | `min_genome_length` | 1700 | Based on the Hepatitis delta (of the _Deltavirus_ genus) genome, the smallest viral genome as of 2024-04-11 (1,700 bp) |
    | `max_genome_length` | 2673870 | Based on the _Pandoravirus salinus_ genome, the largest known viral genome (2,473,870 bp), plus an additional 200 kbp to cater for potential extra genomic material |
    | `min_coverage` | 10 | A bare-minimum average per base coverage across the genome required for genome characterization A Higher coverage would be required for high-quality phylogenetics. |
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
    | Variable  | Rationale |
    | --- | --- | --- |
    | `skip_screen` | false | Set to true to skip the read screen from running. If you set this value to true, please provide a value for the theiaeuk_illumina_pe `genome_length` optional input, OR set the theiaeuk_illumina_pe `call_rasusa` optional input to false. Otherwise RASUSA will attempt to downsample to an expected genome size of 0 bp, and the workflow will fail. |
    | `min_reads` | 3000 | Calculated from the minimum number of base pairs required for 20x coverage of the _Hansenula polymorpha_ genome, the smallest fungal genome as of 2015-04-02 (8.97 Mbp), divided by 300 (the longest Illumina read length) |
    | `min_basepairs` | 45000000 | Should be greater than 10x coverage of _Hansenula polymorpha_, the smallest fungal genome as of 2015-04-02 (8.97 Mbp)  |
    | `min_genome_length` | 9000000 | Based on the _Hansenula polymorpha_  genome - the smallest fungal genome as of 2015-04-02 (8.97 Mbp) |
    | `max_genome_length` | 178000000 | Based on the _Cenococcum geophilum_  genome, the largest pathogenic fungal genome (177.57 Mbp), plus an additional 2 Mbp to cater for potential extra genomic material |
    | `min_coverage` | 10 | A bare-minimum average per base coverage across the genome required for genome characterization. Higher coverage would be required for high-quality phylogenetics.|
    | `min_proportion` | 40 | Neither read1 nor read2 files should have less than 40% of the total number of reads. For paired-end data only. |
<!-- endif -->

    !!! techdetails "Screen Technical Details"
        
        There is a single WDL task for read screening. The `screen` task is run twice, once for raw reads and once for clean reads.
        
        |  | Links |
        | --- | --- |
        | Task | [task_screen.wdl (PE sub-task)](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/comparisons/task_screen.wdl#L3)<br>[task_screen.wdl (SE sub-task)](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/comparisons/task_screen.wdl#L147) |
