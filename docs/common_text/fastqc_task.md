??? task "`FastQC`: Read Quantification"
    `FastQC` quantifies the forward and reverse reads in FASTQ files. For paired-end data, it also provide the total number of read pairs. This task is run once with raw reads as input and once with clean reads as input. If QC has been performed correctly, you should expect **fewer** clean reads than raw reads.

    This tool also provides a graphical visualization of the read quality.
      
    !!! techdetails "`FastQC` Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_fastqc.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_fastqc.wdl) |
        | Software Source Code | [FastQC on Github](https://github.com/s-andrews/FastQC) |
        | Software Documentation | [FastQC Website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) |
