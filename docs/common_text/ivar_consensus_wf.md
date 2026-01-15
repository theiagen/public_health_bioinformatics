
??? task "`ivar_consensus`: Alignment, Consensus, Variant Detection, and Assembly Statistics"
    `iVar Consensus` is a sub-workflow within TheiaCoV that performs reference-based consensus assembly using the [iVar](https://andersen-lab.github.io/ivar/html/index.html) tool by Nathan Grubaugh from the Andersen lab.

{{ include_md("common_text/bwa_task.md", indent=4, condition="theiacov")}}

{{ include_md("common_text/ivar_trim_task.md", indent=4) }}

{{ include_md("common_text/assembly_metrics_task.md", indent=4, condition="theiacov") }}

{{ include_md("common_text/ivar_consensus_task.md", indent=4, condition="theiacov") }}

{{ include_md("common_text/ivar_variants_task.md", indent=4, condition="theiacov") }}

    !!! techdetails "iVar Consensus Technical Details"
        |  | Links |
        | --- | --- |
        | Subworkflow | [wf_ivar_consensus.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_ivar_consensus.wdl) |
