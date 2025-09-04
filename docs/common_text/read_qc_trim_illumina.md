??? task "`read_QC_trim`: Read Quality Trimming, Adapter Removal, Quantification, and Identification"

    `read_QC_trim` is a sub-workflow that removes low-quality reads, low-quality regions of reads, and sequencing adapters to improve data quality. It uses a number of tasks, described below. The differences between the PE and SE versions of the `read_QC_trim` sub-workflow lie in the default parameters, the use of two or one input read file(s), and the different output files.

<!-- if: theiacov|freyja|theiaviral -->
{{ include_md("common_text/ncbi_scrub_task.md", indent=4, replacements={"??? task": "??? toggle"}) }}
<!-- endif -->

    ??? toggle "Read quality trimming"

        Either `trimmomatic` or `fastp` can be used for read-quality trimming. Trimmomatic is used by default. Both tools trim low-quality regions of reads with a sliding window (with a window size of `trim_window_size`), cutting once the average quality within the window falls below `trim_quality_trim_score`. They will both discard the read if it is trimmed below `trim_minlen`. 

        ??? dna "`read_processing` input parameter"
            This input parameter accepts either `trimmomatic` or `fastp` as an input to determine which tool should be used for read quality trimming. This is set to `trimmomatic` by default.

            If the `fastp` option is selected, see below for table of default parameters.

            ??? toggle "`fastp` default read-trimming parameters"
                | **Parameter** | **Explanation** |
                | --- | --- |
                | -g | enables polyG tail trimming |
                | -5 20 | enables read end-trimming |
                | -3 20 | enables read end-trimming |
                | --detect_adapter_for_pe | enables adapter-trimming **only for paired-end reads** |

                Additional arguments can be passed using the `fastp_args` optional parameter.

        !!! techdetails "Trimmomatic and fastp Technical Details"
            |  | Links |
            | --- | --- |
            | Task | [task_trimmomatic.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_trimmomatic.wdl)<br>[task_fastp.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_fastp.wdl) |
            | Software Source Code | [Trimmomatic](https://github.com/usadellab/Trimmomatic)<br>[fastp on Github](https://github.com/OpenGene/fastp) |
            | Software Documentation | [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)<br>[fastp](https://github.com/OpenGene/fastp) |
            | Original Publication(s) | [Trimmomatic: a flexible trimmer for Illumina sequence data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4103590/)<br>[fastp: an ultra-fast all-in-one FASTQ preprocessor](https://academic.oup.com/bioinformatics/article/34/17/i884/5093234?login=false) |

    ??? toggle "Adapter removal"

        The `BBDuk` task removes adapters from sequence reads. To do this:

        - [Repair](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/repair-guide/) from the [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) package reorders reads in paired fastq files to ensure the forward and reverse reads of a pair are in the same position in the two fastq files.
        - [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/)  (*"Bestus Bioinformaticus" Decontamination Using Kmers*) is then used to trim the adapters and filter out all reads that have a 31-mer match to [PhiX](https://emea.illumina.com/products/by-type/sequencing-kits/cluster-gen-sequencing-reagents/phix-control-v3.html), which is commonly added to Illumina sequencing runs to monitor and/or improve overall run quality.

        ??? question "What are adapters and why do they need to be removed?"
            Adapters are manufactured oligonucleotide sequences attached to DNA fragments during the library preparation process. In Illumina sequencing, these adapter sequences are required for attaching reads to flow cells. You can read more about Illumina adapters [here](https://emea.support.illumina.com/bulletins/2020/06/illumina-adapter-portfolio.html). For genome analysis, it's important to remove these sequences since they're not actually from your sample. If you don't remove them, the downstream analysis may be affected.

        !!! techdetails "BBDuk Technical Details"
            |  | Links |
            | --- | --- |
            | Task | [task_bbduk.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_bbduk.wdl) |
            | Software Source Code | [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) |
            | Software Documentation | [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) |
        
    ??? toggle "Read Quantification"

        There are two methods for read quantification to choose from: [`fastq-scan`](https://github.com/rpetit3/fastq-scan) (default) or [`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Both quantify the forward and reverse reads in FASTQ files. For paired-end data, they also provide the total number of read pairs. This task is run once with raw reads as input and once with clean reads as input. If QC has been performed correctly, you should expect **fewer** clean reads than raw reads. `fastqc` also provides a graphical visualization of the read quality.

        ??? dna "`read_qc` input parameter"
            This input parameter accepts either `"fastq_scan"` or `"fastqc"` as an input to determine which tool should be used for read quantification. This is set to `"fastq-scan"` by default.

        !!! techdetails "fastq-scan and FastQC Technical Details"
            |  | Links |
            | --- | --- |
            | Task | [task_fastq_scan.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_fastq_scan.wdl)<br>[task_fastqc.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_fastqc.wdl") |
            | Software Source Code | [fastq-scan on Github](https://github.com/rpetit3/fastq-scan)<br>[fastqc on Github](https://github.com/s-andrews/FastQC) |
            | Software Documentation | [fastq-scan](https://github.com/rpetit3/fastq-scan/blob/master/README.md)<br>[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) |

<!-- if: theiaprok|theiameta -->
[{ include_md("common_text/midas_task.md", indent=4) }]
<!-- endif -->
<!-- if: theiacov -->
{{ include_md("common_text/kraken2_task.md", condition="theiacov", indent=4, replacements={"??? task": "??? toggle"}) }}
<!-- endif -->

<!-- if: theiaviral -->
{{ include_md("common_text/host_decontaminate.md", condition="theiaviral", indent=4, replacements={"??? task": "??? toggle"}) }}

{{ include_md("common_text/kraken2_task.md", condition="theiaviral", indent=4) }}

{{ include_md("common_text/krakentools_task.md", condition="theiaviral", indent=4, replacements={'??? task "`krakentools`"' : '??? toggle "Read Extraction"'}) }}
<!-- endif -->
