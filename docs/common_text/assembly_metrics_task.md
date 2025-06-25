??? task "`assembly_metrics`"

    The `assembly_metrics` task generates mapping statistics from a BAM file. It uses samtools to generate a summary of the mapping statistics, which includes coverage, depth, average base quality, average mapping quality, and other relevant metrics.

    !!! techdetails "`assembly_metrics` Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_assembly_metrics.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_assembly_metrics.wdl) |
        | Software Source Code | [samtools on GitHub](https://github.com/samtools/samtools) |
        | Software Documentation | [samtools](https://www.htslib.org/doc/samtools.html) |
        | Original Publication(s) | [The Sequence Alignment/Map format and SAMtools](https://doi.org/10.1093/bioinformatics/btp352)<br>[Twelve Years of SAMtools and BCFtools](https://doi.org/10.1093/gigascience/giab008) |