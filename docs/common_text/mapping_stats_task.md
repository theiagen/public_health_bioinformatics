---
title: Task Fragment `mapping_stats`
fragment: true
---
??? task "`mapping_stats`: Read Mapping Statistics"
    The Read Mapping Statistics task generates mapping statistics from a BAM file. It uses SAMtools to generate a summary of the mapping statistics, which includes coverage, depth, average base quality, average mapping quality, and other relevant metrics. These statistics are also reported on a per sequence basis.

    !!! techdetails "Read Mapping Statistics Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_mapping_stats.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_mapping_stats.wdl) |
        | Software Source Code | [SAMtools on GitHub](https://github.com/samtools/samtools) |
        | Software Documentation | [SAMtools](https://www.htslib.org/doc/samtools.html) |
        | Original Publication(s) | [The Sequence Alignment/Map format and SAMtools](https://doi.org/10.1093/bioinformatics/btp352)<br>[Twelve Years of SAMtools and BCFtools](https://doi.org/10.1093/gigascience/giab008) |
