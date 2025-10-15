??? task "`kraken_parser`: Parses Kraken Reports"
    
    `kraken_parser` lightens the computation load by taking the input taxon ID list and comparing it to the taxon IDs identified by Kraken2 in the `kraken_report_clean` output file. Only taxon IDs that were found by Kraken are used in the scatter portion of the workflow, which lowers the number of scatter shards the workflow requires.

    !!! techdetails "Find Files Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_kraken_parser.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/file_handling/task_kraken_parser.wdl) |