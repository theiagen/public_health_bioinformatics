??? task "`kraken_parser`: Parses Kraken Reports"
    
    TheiaViral_Panel uses a large input set of taxon IDs for read classification. `kraken_parser` lightens the computational load this creates by taking the input taxon ID list and the clean Kraken report and determines which of the input IDs have been identified by Kraken. A new taxon list is returned to be used in the scatter portion of the workflow lowering the number of scatter shards the workflow spins up. 

    !!! techdetails "Find Files Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_find_files.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/file_handling/task_kraken_parser.wdl) |