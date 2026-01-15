??? task "`artic_consensus`: Alignment, Primer Trimming, Variant Detection, and Consensus"
    This task runs the `Artic minion` command which is a pipeline with a number of stages, [described in detail in the ARTIC documentation](https://artic.readthedocs.io/en/latest/minion/#0). Briefly, these stages are as follows:

    Input reads are aligned to the appropriate reference and only mapped reads are retained. Alignment post-processing occurs, where primers are removed and various trimming steps are undertaken. Variants are detected, and a consensus assembly file is generated.

    Please note that the Medaka model is set by default to `"r941_min_high_g360"` which may not be suitable for your sequencing data. Please be sure to change this parameter if needed.

    !!! info ""
        Read-trimming is performed on raw read data generated on the ClearLabs instrument and thus not a required step in the TheiaCoV_ClearLabs workflow.

    !!! techdetails "Artic Consensus Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_artic_consensus.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_artic_consensus.wdl) |
        | Software Source Code | [ARTIC on GitHub](https://github.com/artic-network/fieldbioinformatics/) |
        | Software Documentation | [ARTIC Documentation](https://artic.readthedocs.io/en/latest/?badge=latest) |
