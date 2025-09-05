??? task "`fastq-scan`: Read Quantification"
    `fastq-scan` quantifies the forward and reverse reads in FASTQ files. For paired-end data, it also provide the total number of read pairs. This task is run once with raw reads as input and once with clean reads as input. If QC has been performed correctly, you should expect **fewer** clean reads than raw reads.

    !!! techdetails "`fastq-scan` Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_fastq_scan.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_fastq_scan.wdl)|
        | Software Source Code | [fastq-scan on GitHub](https://github.com/rpetit3/fastq-scan) |
        | Software Documentation | [fastq-scan on GitHub](https://github.com/rpetit3/fastq-scan/blob/master/README.md) |
