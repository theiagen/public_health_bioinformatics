??? task "`arln_stats`: Quality Assessment for ARLN (optional)"

    To activate this task, set `call_arln_stats` to `true`.

    The `arln_stats` task will provide the user with Antimicrobial Resistance Laboratory Network ([ARLN](https://www.cdc.gov/antimicrobial-resistance-laboratory-networks/php/about/domestic.html))-compliant statistics used for PASS/FAIL assessment, such as assembly ratio, percent GC statistics, and the percent Q30 of raw and cleaned reads.

    The Q30 statistics are calculated via a Python script ([q30.py](https://github.com/theiagen/theiagen_docker_builds/blob/main/arln_stats/1.0.0/q30.py)) that is incorporated within this task's Docker image. The assembly ratio statistic is obtained by parsing an included [NCBI Assembly Stats](https://github.com/CDCgov/phoenix/blob/main/assets/databases/NCBI_Assembly_stats_20240124.txt) text file that is created by the CDC to aggregate NCBI prokaryotic assembly data. The file is sorted by taxon and has statistics for each taxon (see toggle below). This data is used to obtain and calculate the assembly ratio of samples passed through `arln_stats`, much the same as CDC's [Phoenix pipeline](https://github.com/CDCgov/phoenix/wiki/Pipeline-Overview#qc-of-assembled-scaffolds--500bps).

    ??? toggle "NCBI Assembly Stats Explained"
        This aggregated data file is produced by the [CDC](https://github.com/CDCgov/phoenix/wiki/Pipeline-Overview#qc-of-assembled-scaffolds--500bps) and the columns are as follows: 

        - Species
        - Assembly_Size_Min
        - Assembly_Size_Max
        - Assembly_vMedian
        - Assembly_Size_Mean
        - Assembly_Size_StDev
        - Assembly_count
        - GC_Min
        - GC_Max
        - GC_Median
        - GC_Mean
        - GC_Stdev
        - GC_count
        - CDS_Min
        - CDS_Max
        - CDS_Median
        - CDS_Mean
        - CDS_Stdev
        - CDS_count
        - Consensus_TAXID
    
    !!! techdetails "arln_stats Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_arln_stats.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/data_handling/task_arln_stats.wdl) |