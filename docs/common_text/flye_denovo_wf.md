---
title: Workflow Fragment `flye_denovo`
fragment: true
---
??? task "`flye_denovo`: _De novo_ Assembly"
    _De novo_ assembly is the process or product of attempting to reconstruct a genome from scratch (without prior knowledge of the genome).

{{ include_md("common_text/porechop_task.md", indent=4) }}
{{ include_md("common_text/flye_task.md", indent=4) }}
{{ include_md("common_text/bandage_task.md", indent=4) }}
{{ include_md("common_text/polypolish_task.md", indent=8) }}

    !!! dna ""
        To skip polishing,, set `skip_polishing` to `true`. These three modules are mututally exlcusive.

{{ include_md("common_text/medaka_task.md", indent=8) }}
{{ include_md("common_text/dorado_task.md", indent=8) }}
{{ include_md("common_text/racon_task.md", indent=8) }}

{{ include_md("common_text/filter_contigs_task.md", indent=4, condition="flye") }}

{{ include_md("common_text/dnaapler_task.md", indent=4) }}

    !!! techdetails "Flye-Denovo Technical Details"
        |  | Links
        | --- | --- |
        | Subworkflow | [wf_flye_denovo.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_flye_denovo.wdl) |
