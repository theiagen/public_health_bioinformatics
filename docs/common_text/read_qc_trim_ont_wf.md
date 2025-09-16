??? task "`read_QC_trim_ont`: Read Quality Trimming, Quantification, and Identification"
    `read_QC_trim_ont` is a sub-workflow that filters low-quality reads and trims low-quality regions of reads. It uses several tasks, described below.

<!-- if: theiacov|freyja -->
{{ include_md("common_text/ncbi_scrub_task.md", indent=4) }}

{{ include_md("common_text/artic_guppyplex_task.md", indent=4) }}

{{ include_md("common_text/kraken2_task.md", indent=4, condition="theiacov") }}
  
<!-- endif -->
  
<!-- if: theiaprok -->
    !!! dna "A note on estimated genome length"

        By default, the estimated genome length is set to 5 Mb, which is around 0.7 Mb higher than the average bacterial genome length, according to [the information of thousands of NCBI bacterial assemblies collated here](https://github.com/CDCgov/phoenix/blob/717d19c19338373fc0f89eba30757fe5cfb3e18a/assets/databases/NCBI_Assembly_stats_20240124.txt). This estimate can be overwritten by the user and is used by `Rasusa`.
<!-- endif -->

<!-- if: theiaprok -->

{{ include_md("common_text/rasusa_task.md", indent=4, condition="ont") }}

{{ include_md("common_text/nanoq_task.md", indent=4) }}

{{ include_md("common_text/kraken2_task.md", indent=4, condition="ont") }}
<!-- endif -->

{{ include_md("common_text/nanoplot_task.md", indent=4, condition="ont") }}

    !!! techdetails "read_QC_trim_ont Technical Details"
        |  | Links |
        | --- | --- |
        | Subworkflow | [wf_read_QC_trim_ont.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_read_QC_trim_ont.wdl) |
