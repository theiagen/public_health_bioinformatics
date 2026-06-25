---
title: Task Fragment `contaminant_check`
fragment: true
---
??? task "`Contaminant_Check`: Contaminant Sequence Status Check"
    The Contaminant Check task outputs a pass/fail status based on if contaminant/host sequences pass thresholds for breadth of coverage, depth of coverage, and number of reads mapped. This task is activated by inputting a comma-delimited string of `expected_sequences`, which match sequence headers in the inputted contaminant/host FASTA. Alternatively, a JSON, `expected_sequences_json`, can be inputted, which relates contaminant sequence set names to a list of contaminant sequences; in this instance, `expected_sequences` is populated with the set name. `expected_sequences_json` must be formatted as follows:

    ```json
    {"<SET_NAME1>": ["<CONTAMINANT_SEQUENCE_1>", "<CONTAMINANT_SEQUENCE_n">]}
    ```
    
    The composite status, `contaminant_check_status`, will report `"PASS"` if expected and unexpected sequences are identified within the `min_expected_seq` and `max_unexpected_seq` thresholds; if not, `"FAIL ..."` is reported depicting which `expected_sequences` failed and why, along with which unexpected sequences were identified. 

    Additionally, the coverage, depth, and number of reads mapped are reported in JSON mappings for the sets of expected and unexpected sequences. 

    !!! techdetails "Contaminant Check Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_contaminant_check.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/contamination/task_contaminant_check.wdl) |
