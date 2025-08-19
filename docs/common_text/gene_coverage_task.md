??? task "`gene_coverage`"

    This task calculates the percent of a region (typically genes) covered above a minimum depth using `samtools` and basic arithmetic. By default, this task runs for SARS-CoV-2 and Mpox, but if a BED file is provided with regions of interest, this task can run for other organisms as well.

    !!! techdetails "Gene Coverage Technical Details"        
        |  | Links |
        | --- | --- |
        | Task | [task_gene_coverage.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_gene_coverage.wdl) |
        | Software Source Code | [SAMtools on GitHub](https://github.com/samtools/samtools) |
        | Software Documentation | [SAMTools Manual](https://www.htslib.org/doc/samtools.html) |
        | Original Publication(s) | [Twelve years of SAMtools and BCFtools](https://doi.org/10.1093/gigascience/giab008) |
