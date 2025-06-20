??? task "`ARLN_Stats`: Stats required by ARLN such as assembly ratio and Q30% scores (optional)"

    The `arln_stats` task when enabled, will provide the user with not already provided ARLN compliant statistics used for PASS/FAIL assesment such as assembly ratio, GC Percent ratio, assembly and GC zscore, and will calculate the Q30% scores of the raw and cleaned reads, provided that these are utilized by the run workflow. TheiaProk_FASTA will not provide these results, only assembly ratio, and likewise TheiaProk_Illumina_SE and TheiaProk_ONT will only provide R1 Q30% metrics.
    
    The Q30% statistics are calculated via a Python script, [q30.py](https://github.com/theiagen/theiagen_docker_builds/blob/main/arln_stats/1.0.0/q30.py), that is incorporated within the Docker image. The assembly ratio statistic is obtained by parsing an included [NCBI_Stats](../../assets/files/NCBI_Assembly_stats_20240124.txt) text file that is created by the CDC to aggregate NCBI prokayotic assembly data. The file is sorted by taxon and has statistics such as minimum, maximum, and mean assembly length, the assembly ratio, and the standard deviation for each taxon. This data is used to obtain and calculate the assembly ratio of samples passed through ARLN_Stats, much the same as CDC's Phoenix. The NCBI_Stats file will be kept up to date when there are updates to the parent NCBI dataset.
    
    ??? toggle "NCBI_Assembly_stats_20240124.txt Explained"
    
        This aggregated data file is produced by the [CDC](https://github.com/CDCgov/phoenix/wiki/Pipeline-Overview#qc-of-assembled-scaffolds--500bps) and the columns are as follows: 

        - Species
        - Assembly_Size_Min
        - Assembly_ Size_Max
        - Assembly_vMedian
        - Assembly_ Size_Mean
        - Assembly_ Size_StDev
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
    
    !!! techdetails "ARLN_Stats Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_arln_stats.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/data_handling/task_arln_stats.wdl) |