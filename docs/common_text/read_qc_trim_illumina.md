??? task "`read_QC_trim`: Read Quality Trimming, Adapter Removal, Quantification, and Identification"

    `read_QC_trim` is a sub-workflow that removes low-quality reads, low-quality regions of reads, and sequencing adapters to improve data quality. It uses a number of tasks, described below. The differences between the PE and SE versions of the `read_QC_trim` sub-workflow lie in the default parameters, the use of two or one input read file(s), and the different output files.

<!-- if: theiacov|freyja|theiaviral -->
{{ include_md("common_text/ncbi_scrub_task.md", indent=4, replacements={"??? task": "??? toggle"}) }}
<!-- endif -->

    ??? toggle "Read quality trimming"

        !!! dna "`read_processing` with `"trimmomatic"` (default) or `"fastp"`"
            Either `trimmomatic` or `fastp` can be used for read-quality trimming. Trimmomatic is used by default. 
            
            To activate `fastp`, set the `read_processing` input parameter to `"fastp"`. 
            
            These tasks are mutually exclusive.

{{ include_md("common_text/trimmomatic_task.md", indent=8) }}  

{{ include_md("common_text/fastp_task.md", indent=8) }}

    ??? toggle "Adapter removal"

{{ include_md("common_text/bbduk_task.md", indent=8) }}

    ??? toggle "Read Quantification"

        !!! dna "`read_qc` with `"fastq-scan"` (default) or `"fastqc"`"
            Either `fastq-scan` or `fastqc` can be used for read quantification. `fastq-scan` is used by default. 
            
            To activate `fastqc`, set the `read_qc` input parameter to `"fastqc"`.
            
            These tasks are mutually exclusive.

{{ include_md("common_text/fastq_scan_task.md", indent=8) }}

{{ include_md("common_text/fastqc_task.md", indent=8) }}

<!-- if: theiaprok|theiameta -->
{{ include_md("common_text/midas_task.md", indent=4, replacements={"??? task": "??? toggle"}) }}
<!-- endif -->

<!-- if: theiacov -->
{{ include_md("common_text/kraken2_task.md", condition="theiacov", indent=4, replacements={"??? task": "??? toggle"}) }}
<!-- endif -->

<!-- if: theiaviral -->
{{ include_md("common_text/host_decontaminate.md", condition="theiaviral", indent=4, replacements={"??? task": "??? toggle"}) }}

{{ include_md("common_text/kraken2_task.md", condition="theiaviral", indent=4) }}

{{ include_md("common_text/krakentools_task.md", condition="theiaviral", indent=4, replacements={'??? task "`krakentools`"' : '??? toggle "Read Extraction"'}) }}
<!-- endif -->
