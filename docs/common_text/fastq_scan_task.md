??? task "`fastq-scan`"
    
    `fastq-scan` provides read quantification of both the forward and reverse reads in FASTQ files. When reads are paired end `fastq-scan` also provides total number of reads.

<!-- if: theiaviral_panel -->
    When used as part of the `TheiaViral_Panel` workflow, `fastq-scan` is used to quantify the number of reads to a given threshold in order to ensure there is enough extracted data to proceed to assembly.
<!-- endif -->

    !!! techdetails "fastq-scan and FastQC Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_fastq_scan.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_fastq_scan.wdl) |
        | Software Source Code | [fastq-scan on Github](https://github.com/rpetit3/fastq-scan) |
        | Software Documentation | [fastq-scan](https://github.com/rpetit3/fastq-scan/blob/master/README.md) |