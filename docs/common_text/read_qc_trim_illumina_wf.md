??? task "`read_QC_trim`: Read Quality Trimming, Adapter Removal, Quantification, and Identification"
    `read_QC_trim` is a sub-workflow that removes low-quality reads, low-quality regions of reads, and sequencing adapters to improve data quality. It uses a number of tasks, described below. The differences between the PE and SE versions of the `read_QC_trim` sub-workflow lie in the default parameters, the use of two or one input read file(s), and the different output files.

<!-- if: theiacov|freyja -->
{{ include_md("common_text/ncbi_scrub_task.md", indent=4) }}
<!-- endif -->

    !!! dna ""
        By default, `read_processing` is set to `"trimmomatic"`. To use `fastp` instead, set `read_processing` to `"fastp"`. These tasks are mutually exclusive.

{{ include_md("common_text/trimmomatic_task.md", indent=8) }}  

{{ include_md("common_text/fastp_task.md", condition="read_qc_trim", indent=8) }}

{{ include_md("common_text/bbduk_task.md", indent=4) }}

    !!! dna ""
        By default, `read_qc` is set to `"fastq_scan"`. To use `fastqc` instead, set `read_qc` to `"fastqc"`. These tasks are mutually exclusive.

{{ include_md("common_text/fastq_scan_task.md", condition="notclearlabs", indent=8) }}

{{ include_md("common_text/fastqc_task.md", indent=8) }}

{{ include_md("common_text/host_decontaminate_wf.md", indent=4) }}

<!-- if: theiaprok|theiameta -->
{{ include_md("common_text/midas_task.md", indent=4) }}
<!-- endif -->

<!-- if: theiaprok|theiaeuk -->
{{ include_md("common_text/kraken2_task.md", condition="theiaprokillumina", indent=4) }}
<!-- endif -->

<!-- if: theiacov -->
{{ include_md("common_text/kraken2_task.md", condition="theiacov", indent=4) }}
<!-- endif -->

    !!! techdetails "read_QC_trim Technical Details"
        |  | Links |
        | --- | --- |
        | Subworkflow | [wf_read_QC_trim_pe.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_read_QC_trim_pe.wdl)<br>[wf_read_QC_trim_se.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_read_QC_trim_se.wdl) |
