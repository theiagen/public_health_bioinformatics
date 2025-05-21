
??? task "`qc_check`: Check QC Metrics Against User-Defined Thresholds (optional)"

    The `qc_check` task compares generated QC metrics against user-defined thresholds for each metric. This task will run if the user provides a `qc_check_table` TSV file. If all QC metrics meet the threshold, the `qc_check` output variable will read `QC_PASS`. Otherwise, the output will read `QC_NA` if the task could not proceed or `QC_ALERT` followed by a string indicating what metric failed.

<!-- if: theiacov -->
    The `qc_check` task applies quality thresholds according to the specified organism, which should match the _standardized_ `organism` input in the TheiaCoV workflows.
<!-- endif -->
<!-- if: theiaprok|theiaeuk -->
    The `qc_check` task applies quality thresholds according to the sample taxa. The sample taxa is taken from the `gambit_predicted_taxon` value inferred by the GAMBIT module OR can be manually provided by the user using the `expected_taxon` workflow input.
<!-- endif -->

    ??? toggle "Formatting the _qc_check_table.tsv_"
        - The first column of the qc_check_table lists the `organism` that the task will assess and the header of this column must be "**taxon**".
<!-- if: theiaprok|theiaeuk -->
        - Any genus or species can be included as a row of the qc_check_table. However, these taxa must **uniquely** match the sample taxa, meaning that the file can include multiple species from the same genus (Vibrio_cholerae and Vibrio_vulnificus), but not both a genus row and species within that genus (Vibrio and Vibrio cholerae). **The taxa should be formatted with the first letter capitalized and underscores in lieu of spaces.**
<!-- endif -->
        - Each subsequent column indicates a QC metric and lists a threshold for each organism that will be checked. **The column names must exactly match expected values, so we highly recommend copy and pasting the header from the template file below as a starting place.**
    
    ??? toggle "Template _qc_check_table.tsv_ files"    
<!-- if: theiacov -->
        - TheiaCoV_Illumina_PE: [TheiaCoV_Illumina_PE_qc_check_template.tsv](../../assets/files/TheiaCoV_Illumina_PE_qc_check_template.tsv)
<!-- endif -->
<!-- if: theiaprok -->
        - TheiaProk_Illumina_PE: [theiaprok_illumina_pe_qc_check_template.tsv](../../assets/files/TheiaProk_Illumina_PE_qc_check_template.tsv)
        - TheiaProk_FASTA: [theiaprok_fasta_qc_check_template.tsv](../../assets/files/TheiaProk_FASTA_qc_check_template.tsv)
<!-- endif -->
<!-- if: theiaeuk -->
        - TheiaEuk_Illumina_PE_PHB: [theiaeuk_qc_check_template.tsv](../../assets/files/TheiaEuk_qc_check_template.tsv)
<!-- endif -->
        !!! warning "Example Purposes Only"
            The QC threshold values shown in the file above are for example purposes only and should not be presumed to be sufficient for every dataset.

    !!! techdetails "`qc_check` Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_qc_check_phb.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/comparisons/task_qc_check_phb.wdl) |
