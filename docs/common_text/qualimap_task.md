??? task "`QualiMap`: BAM File Quality Assessment"

    QualiMap evaluates the quality of alignment data in BAM files by computing various metrics including coverage distribution, mapping quality, GC content, and various metrics analyzed across the reference. It provides comprehensive quality control reports for next-generation sequencing alignment data.

    This task generates both standard QualiMap reports and custom interactive HTML visualizations for genome coverage and mapping quality across the reference sequence. The results are bundled into a compressed archive for easy download and review, especially since for the QualiMap report to render the pngs correctly, it needs to preserve directory structure. 

    !!! techdetails "QualiMap Technical Details"

        |  | Links |
        | --- | --- |
        | Task | [task_qualimap.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_qualimap.wdl) |
        | Software Source Code | [QualiMap on Bitbucket](https://bitbucket.org/kokonech/qualimap/src/master/) |
        | Software Documentation | [QualiMap Documentation](http://qualimap.conesalab.org/) |
        | Original Publication(s) | [QualiMap: evaluating next-generation sequencing alignment data](https://academic.oup.com/bioinformatics/article/28/20/2678/206551) |
