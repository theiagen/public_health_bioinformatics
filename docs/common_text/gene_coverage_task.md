??? task "`gene_coverage`"

    This task calculates the percent of the gene covered above a minimum depth. By default, it runs for SARS-CoV-2 and MPXV, but if a bed file is provided with regions of interest, this task will be run for other organisms as well.

    !!! techdetails "Gene Coverage Technical Details"        
        |  | Links |
        | --- | --- |
        | Task | [task_gene_coverage.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_gene_coverage.wdl) |
