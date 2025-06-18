??? task "`nanoplot`"

    Nanoplot is used for the determination of mean quality scores, read lengths, and number of reads. This task is run once with raw reads as input and once with clean reads as input. If QC has been performed correctly, you should expect **fewer** clean reads than raw reads.

    !!! techdetails "Nanoplot Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_nanoplot.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_nanoplot.wdl) |
        | Software Source Code | [NanoPlot](https://github.com/wdecoster/NanoPlot) |
        | Software Documentation | [NanoPlot Documentation](https://github.com/wdecoster/NanoPlot/blob/master/README.md) |
        | Original Publication(s) | [NanoPack2: population-scale evaluation of long-read sequencing data](https://academic.oup.com/bioinformatics/article/39/5/btad311/7160911) |