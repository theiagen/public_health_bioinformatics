??? task "`NanoPlot`: Read Quantification"

    NanoPlot is used for the determination of mean quality scores, read lengths, and number of reads. This task is run once with raw reads as input and once with clean reads as input. If QC has been performed correctly, you should expect **fewer** clean reads than raw reads.

<!-- if: ont -->
    While this task currently is run _outside_ of the `read_QC_trim_ont` workflow, it is being included here as it calculates statistics on the read data. This is done so that the actual assembly genome lengths can be used (if an estimated genome length is not provided by the user) to ensure the estimated coverage statistics are accurate.
<!-- endif -->

    !!! techdetails "NanoPlot Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_nanoplot.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_nanoplot.wdl) |
        | Software Source Code | [NanoPlot on GitHub](https://github.com/wdecoster/NanoPlot) |
        | Software Documentation | [NanoPlot Documentation](https://github.com/wdecoster/NanoPlot/blob/master/README.md) |
        | Original Publication(s) | [NanoPack2: population-scale evaluation of long-read sequencing data](https://academic.oup.com/bioinformatics/article/39/5/btad311/7160911) |