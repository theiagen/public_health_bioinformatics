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
  
<!-- if: theiaprok -->
    !!! dna "A note on estimated genome length and Rasusa"

        Previously, `genome_length` was defaulted to 5 Mb, which caused Rasusa to run on every sample. That default has been removed, so Rasusa is now skipped unless the user explicitly requests downsampling by setting one of the following: `genome_length`, `rasusa_num_bases`, `rasusa_fraction_of_reads`, or `rasusa_num_reads`. To downsample by coverage, please provide the genome length of your expected taxa. Downsampling is beneficial when coverage exceeds ~150x as this can increase assembly time, computational costs, and can introduce errors. 

{{ include_md("common_text/rasusa_task.md", indent=4, condition="ont") }}

{{ include_md("common_text/nanoq_task.md", indent=4) }}

{{ include_md("common_text/metabuli_task.md", indent=4, condition="theiaprok") }}
<!-- endif -->

!!! dna "Nanoplot and genome length"

    If `genome_length` is not provided Nanoplot will utilize the Quast assembly length to calculate the estimated coverage. 

{{ include_md("common_text/nanoplot_task.md", indent=4, condition="ont") }}

    !!! techdetails "read_QC_trim_ont Technical Details"
        |  | Links |
        | --- | --- |
        | Subworkflow | [wf_read_QC_trim_ont.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_read_QC_trim_ont.wdl) |
