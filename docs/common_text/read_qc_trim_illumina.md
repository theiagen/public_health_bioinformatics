??? task "`read_QC_trim`: Read Quality Trimming, Adapter Removal, Quantification, and Identification"

    `read_QC_trim` is a sub-workflow that removes low-quality reads, low-quality regions of reads, and sequencing adapters to improve data quality. It uses a number of tasks, described below. The differences between the PE and SE versions of the `read_QC_trim` sub-workflow lie in the default parameters, the use of two or one input read file(s), and the different output files.

<!-- if: theiacov|freyja|theiaviral -->
{{ include_md("common_text/ncbi_scrub_task.md", indent=4, replacements={"??? task": "??? toggle"}) }}
<!-- endif -->

    ??? toggle "Read quality trimming"

        Either `trimmomatic` or `fastp` can be used for read-quality trimming. Trimmomatic is used by default. Both tools trim low-quality regions of reads with a sliding window (with a window size of `trim_window_size`), cutting once the average quality within the window falls below `trim_quality_trim_score`. They will both discard the read if it is trimmed below `trim_minlen`. 

        ??? dna "`read_processing`"
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

        ??? dna "`read_qc`"
            This input parameter accepts either `"fastq_scan"` or `"fastqc"` as an input to determine which tool should be used for read quantification. This is set to `"fastq-scan"` by default.

        !!! techdetails "fastq-scan and FastQC Technical Details"
            |  | Links |
            | --- | --- |
            | Task | [task_fastq_scan.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_fastq_scan.wdl)<br>[task_fastqc.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_fastqc.wdl") |
            | Software Source Code | [fastq-scan on Github](https://github.com/rpetit3/fastq-scan)<br>[fastqc on Github](https://github.com/s-andrews/FastQC) |
            | Software Documentation | [fastq-scan](https://github.com/rpetit3/fastq-scan/blob/master/README.md)<br>[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) |

<!-- if: theiaprok|theiameta -->
    ??? toggle "Read Identification with MIDAS (optional)"
        The `MIDAS` task is for the identification of reads to detect contamination with non-target taxa. This task is optional and turned off by default. It can be used by setting the `call_midas` input variable to `true`.

        The MIDAS tool was originally designed for metagenomic sequencing data but has been co-opted for use with bacterial isolate WGS methods. It can be used to detect contamination present in raw sequencing data by estimating bacterial species abundance in bacterial isolate WGS data. If a secondary genus is detected above a relative frequency of 0.01 (1%), then the sample should fail QC and be investigated further for potential contamination.

        This task is similar to those used in commercial software, BioNumerics, for estimating secondary species abundance.

        ??? toggle "How are the MIDAS output columns determined?"
            
            Example MIDAS report in the `midas_report` column:
            
            | species_id | count_reads | coverage | relative_abundance |
            | --- | --- | --- | --- |
            | Salmonella_enterica_58156 | 3309 | 89.88006645 | 0.855888033 |
            | Salmonella_enterica_58266 | 501 | 11.60606061 | 0.110519371 |
            | Salmonella_enterica_53987 | 99 | 2.232896237 | 0.021262881 |
            | Citrobacter_youngae_61659 | 46 | 0.995216227 | 0.009477003 |
            | Escherichia_coli_58110 | 5 | 0.123668877 | 0.001177644 |
            
            MIDAS report column descriptions:
            
            - species_id: species identifier
            - count_reads: number of reads mapped to marker genes
            - coverage: estimated genome-coverage (i.e. read-depth) of species in metagenome
            - relative_abundance: estimated relative abundance of species in metagenome
            
            The value in the `midas_primary_genus` column is derived by ordering the rows in order of "relative_abundance" and identifying the genus of top species in the "species_id" column (Salmonella). The value in the `midas_secondary_genus` column is derived from the genus of the second-most prevalent genus in the "species_id" column (Citrobacter). The `midas_secondary_genus_abundance` column is the "relative_abundance" of the second-most prevalent genus (0.009477003). The `midas_secondary_genus_coverage` is the "coverage" of the second-most prevalent genus (0.995216227).

        **MIDAS Reference Database Overview**

        The **MIDAS reference database** is a comprehensive tool for identifying bacterial species in metagenomic and bacterial isolate WGS data. It includes several layers of genomic data, helping detect species abundance and potential contaminants.

        !!! dna "Key Components of the MIDAS Database"
            1. **Species Groups**: 
                - MIDAS clusters bacterial genomes based on 96.5% sequence identity, forming over 5,950 species groups from 31,007 genomes. These groups align with the gold-standard species definition (95% ANI), ensuring highly accurate species identification.

            2. **Genomic Data Structure**:
                - _Marker Genes_: Contains 15 universal single-copy genes used to estimate species abundance.
                - _Representative Genome_: Each species group has a selected representative genome, which minimizes genetic variation and aids in accurate SNP identification.
                - _Pan-genome_: The database includes clusters of non-redundant genes, with options for multi-level clustering (e.g., 99%, 95%, 90% identity), enabling MIDAS to identify gene content within strains at various clustering thresholds.

            3. **Taxonomic Annotation**: 
                - Genomes are annotated based on consensus Latin names. Discrepancies in name assignments may occur due to factors like unclassified genomes or genus-level ambiguities.

        ---

        **Using the Default MIDAS Database**

        TheiaProk and TheiaEuk use the pre-loaded MIDAS database in Terra (see input table for current version) by default for bacterial species detection in metagenomic data, requiring no additional setup.

        !!! tip "Create a Custom MIDAS Database"
            Users can also build their own custom MIDAS database if they want to include specific genomes or configurations. This custom database can replace the default MIDAS database used in Terra. To build a custom MIDAS database, follow the [MIDAS GitHub guide on building a custom database](https://github.com/snayfach/MIDAS/blob/master/docs/build_db.md). Once the database is built, users can upload it to a Google Cloud Storage bucket or Terra workkspace and provide the link to the database in the `midas_db` input variable.

        ---

        !!! techdetails "MIDAS Technical Details"
            |  | Links |
            | --- | --- |
            | Task | [task_midas.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/contamination/task_midas.wdl) |
            | Software Source Code | [MIDAS](https://github.com/snayfach/MIDAS)
            | Software Documentation | [MIDAS](https://github.com/snayfach/MIDAS)
            | Original Publication(s) | [An integrated metagenomics pipeline for strain profiling reveals novel patterns of bacterial transmission and biogeography](https://pubmed.ncbi.nlm.nih.gov/27803195) |

<!-- endif -->
<!-- if: theiacov -->
{{ include_md("common_text/kraken2_task.md", condition="theiacov", indent=4, replacements={"??? task": "??? toggle"}) }}
<!-- endif -->

<!-- if: theiaviral -->
{{ include_md("common_text/host_decontaminate.md", condition="theiaviral", indent=4, replacements={"??? task": "??? toggle"}) }}

{{ include_md("common_text/kraken2_task.md", condition="theiaviral", indent=4) }}

{{ include_md("common_text/krakentools_task.md", condition="theiaviral", indent=4, replacements={'??? task "`krakentools`"' : '??? toggle "Read Extraction"'}) }}
<!-- endif -->
