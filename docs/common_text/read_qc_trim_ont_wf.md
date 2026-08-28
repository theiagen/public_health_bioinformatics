---
title: Workflow Fragment `read_qc_trim_ont`
fragment: true
---
??? task "`read_QC_trim_ont`: Read Quality Trimming, Quantification, and Identification"
    `read_QC_trim_ont` is a sub-workflow that filters low-quality reads and trims low-quality regions of reads. It uses several tasks, described below.

<!-- if: theiacov|freyja -->
{{ include_md("common_text/ncbi_scrub_task.md", indent=4) }}

{{ include_md("common_text/artic_guppyplex_task.md", indent=4) }}

{{ include_md("common_text/metabuli_task.md", indent=4, condition="theiacov") }}
<!-- endif -->

<!-- if: theiaprok|theiaeuk -->
{{ include_md("common_text/rasusa_task.md", indent=4, condition="ont") }}
{{ include_md("common_text/nanoq_task.md", indent=4) }}
{{ include_md("common_text/metabuli_task.md", indent=4, condition="theiaprok") }}
<!-- endif -->

{{ include_md("common_text/nanoplot_task.md", indent=4, condition="ont") }}

    !!! techdetails "read_QC_trim_ont Technical Details"
        |  | Links |
        | --- | --- |
        | Subworkflow | [wf_read_QC_trim_ont.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_read_QC_trim_ont.wdl) |
