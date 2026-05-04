
??? task "`contaminant_check`: Contaminant Sequence Identification Status Check"

    `Contaminant Check` determines if contaminant/host sequences pass thresholds for breadth of coverage, depth of coverage, and number of reads mapped, and outputs an overall pass/fail status based on the composite results. This task is activated by inputting a comma-delimited string of `expected_sequences`, which match sequence headers in the inputted contaminant/host FASTA. Each sequence from the previously inputted contaminant/host FASTA is checked for sufficient read mapping statistics. 
    
    The composite status, `contaminant_check_status`, will report `"PASS"` if expected and unexpected sequences are identified within the `min_unexpected_seq` and `max_unexpected_seq` thresholds; if not, `"FAIL ..."` is reported depicting which `expected_sequences` failed and why, and which unexpected sequences were identified. 

    Additionally, the coverage, depth, and number of reads mapped are reported in JSON mappings for the sets of expected and unexpected sequences. 

    !!! techdetails "contaminant_check Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_contaminant_check.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/contamination/task_contaminant_check.wdl) |
