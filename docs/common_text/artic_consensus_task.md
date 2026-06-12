---
title: Task Fragment `artic_consensus`
fragment: true
---
??? task "`artic consensus`: Alignment, Primer Trimming, Variant Detection, and Consensus"
    This task runs the `artic minion` command which is a pipeline with a number of stages, [described in detail in the ARTIC documentation](https://artic.readthedocs.io/en/latest/minion/#0). Briefly, these stages are as follows:

    Input reads are aligned to the appropriate reference and only mapped reads are retained. Alignment post-processing occurs, where primers are removed and various trimming steps are undertaken. Variants are detected, and a consensus assembly file is generated.

    Please note that the Clair3 model is set by default to `"r1041_e82_400bps_sup_v500"` which may not be suitable for your sequencing data. Please be sure to change this parameter if needed.

    !!! info "Primer BED formatting"
        ARTIC has stringent primer BED formatting standards that are depicted [here](https://chrisgkent.github.io/primalbedtools/). The "sequence" and "primerAttributes" fields are NOT necessary for this task, though a comprehensive file is depicted here:
        
        ```
        # chrom     start   end     primername              pool    strand  sequence            primerAttributes       
        MN908947.3  47      78      SARS-CoV-2_1_LEFT_1     1       +       CTCTTGTAGATCTT...   pw=1.0;ps=100
        ```
        
        Some other noteworthy formatting standards include:
        
        - alternative primer schemes must contain an *underscored* suffix, e.g. "_alt#"
        - versions must be *hyphenated* because periods are not permitted e.g. "HIV-v2-0"
        - primer pools *must* be consistent with the primer name's number assignment, e.g. "SARS-CoV-2_1_LEFT" must be in the same primer pool as "SARS-CoV-2_1_RIGHT"

    !!! info "ClearLabs"
        Read-trimming is performed on raw read data generated on the ClearLabs instrument and thus not a required step in the TheiaCoV_ClearLabs workflow.

    !!! techdetails "Artic Consensus Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_artic_consensus.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_artic_consensus.wdl) |
        | Software Source Code | [ARTIC on GitHub](https://github.com/artic-network/fieldbioinformatics/) |
        | Software Documentation | [ARTIC Documentation](https://artic.readthedocs.io/en/latest/?badge=latest) |
