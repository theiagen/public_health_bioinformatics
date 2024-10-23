# TheiaMeta Panel

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Genomice Characterization](../../workflows_overview/workflows_type.md/#genomic_characterization) | [Viral](../../workflows_overview/workflows_kingdom.md/#viral) | PHB v2.X.X | Yes | Sample-level |

## TheiaMeta_Panel_Illumina_PE_PHB

TheiaMeta_Panel was created initially for the Illumina Viral Surveillance Panel; however, it can be used for any panel that is sequenced using Illumina paired-end reads if the appropriate taxon IDs are provided. TheiaMeta_Panel performs taxonomic binning, and then assembles the bins into contigs. If the contigs are associated with a supported organism, genomic characterization will be performed.

### Inputs

<div class="searchable-table" markdown="1">

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| theiameta_panel_illumina_pe | **read1** | File | The forward Illumina read in FASTQ file format (compression optional)  | | Required |
| theiameta_panel_illumina_pe | **read2** | File | The reverse Illumina read in FASTQ file format (compression optional) | | Required |
| theiameta_panel_illumina_pe | **samplename** | String | The name of the sample being analyzed | | Required |
| theiameta_panel_illumina_pe | **taxon_ids** | Array[Int] | The taxon IDs to be used for taxonomic binning | | Required |
| fastq_scan_binned | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional |
| fastq_scan_binned | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| fastq_scan_binned | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1 | Optional |
| fastq_scan_binned | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional |
| gather_scatter | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| gather_scatter | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 50 | Optional |
| gather_scatter | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16 | Optional |
| gather_scatter | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| kraken2 | **classified_out** | String | Allows user to rename the classified FASTQ files output. Must include .fastq as the suffix | classified#.fastq | Optional |
 kraken2 | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional |
| kraken2 | **disk_size** | Int | GB of storage to request for VM used to run the kraken2 task. Increase this when using large (>30GB kraken2 databases such as the "k2_standard" database) | 100 | Optional |
| kraken2 | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/kraken2:2.1.2-no-db | Optional |
| kraken2 | **kraken2_args** | String | Allows a user to supply additional kraken2 command-line arguments |  | Optional |
| kraken2 | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 32 | Optional |
| kraken2 | **unclassified_out** | String | Allows user to rename unclassified FASTQ files output. Must include .fastq as the suffix | unclassified#.fastq | Optional |
| krakentools | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional |
| krakentools | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| krakentools | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/krakentools:d4a2fbe| Optional |
| krakentools | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 4 | Optional |
| metaspades | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional |
| metaspades | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| metaspades | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/metaspades:3.15.3 | Optional |
| metaspades | **kmers** | String | The k-mer list to use; if not provided, the value is automatically set | | Optional |
| metaspades | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 32 | Optional |
| metaspades | **metaspades_opts** | String | Additional arguments to pass on to the metaspades command | | Optional |
| metaspades | **phred_offset** | Int | The PHRED quality offset of the input reads; can be either 33 or 64 | 33 | Optional |
| minimap2_assembly_correction | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| minimap2_assembly_correction | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| minimap2_assembly_correction | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/minimap2:2.22 | Optional |
| minimap2_assembly_correction | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| pilon | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| pilon | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| pilon | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/pilon:1.24 | Optional |
| pilon | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| quast | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| quast | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| quast | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/quast:5.0.2 | Optional |
| quast | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| read_QC_trim | **adapters** | File | A file containing the sequence of the adapters used during library preparation, used in the BBDuk task |  | Optional |
| read_QC_trim | **bbduk_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| read_QC_trim | **call_kraken** | Boolean | Set to true to launch Kraken2; if true, you must provide a kraken_db | FALSE | Optional |
| read_QC_trim | **call_midas** | Boolean | Set to true to launch Midas | TRUE | Optional |
| read_QC_trim | **fastp_args** | String | Additional arguments to pass to fastp | "--detect_adapter_for_pe -g -5 20 -3 20 | Optional |
| read_QC_trim | **midas_db** | File | Midas database file | gs://theiagen-large-public-files-rp/terra/theiaprok-files/midas/midas_db_v1.2.tar.gz | Optional |
| read_QC_trim | **phix** | File | A file containing the phix used during Illumina sequencing; used in the BBDuk task |  | Optional |
| read_QC_trim | **read_processing** | String | Read trimming software to use, either "trimmomatic" or "fastp" | trimmomatic | Optional |
| read_QC_trim | **read_qc** | String | Allows the user to decide between fastq_scan (default) and fastqc for the evaluation of read quality. | fastq_scan | Optional |
| read_QC_trim | **trim_min_length** | Int | The minimum length of each read after trimming | 75 | Optional | 
| read_QC_trim | **trim_primers** | Boolean | A True/False option that determines if primers should be trimmed. | TRUE | Optional |
| read_QC_trim | **trim_quality_min_score** | Int | The minimum quality score to keep during trimming | 30 | Optional |
| read_QC_trim | **trim_window_size** | Int | Specifies window size for trimming (the number of bases to average the quality across) | 4 | Optional |
| read_QC_trim | **trimmomatic_args** | String | Additional arguments to pass to trimmomatic | -phred33 | Optional |
| sort_bam_assembly_correction | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| sort_bam_assembly_correction | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| sort_bam_assembly_correction | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/samtools:1.17 | Optional |
| sort_bam_assembly_correction | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| theiameta_panel_illumina_pe | **kraken2_db** | File | A Kraken2 database in .tar.gz format | gs://theiagen-large-public-files-rp/terra/databases/kraken2/k2_viral_20240112.tar.gz | Optional | 
| theiameta_panel_illumina_pe | **minimum_read_number** | Int | The minimum number of reads in order to attempt assembly on a bin of reads | 1000 | Optional |
| version_capture | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) | | Optional |

</div>

### Workflow Tasks

??? task "`read_QC_trim`: Read Quality Trimming, Adapter Removal, Quantification, and Identification"

    `read_QC_trim` is a sub-workflow within TheiaMeta that removes low-quality reads, low-quality regions of reads, and sequencing adapters to improve data quality. It uses a number of tasks, described below.

    **Read quality trimming**

    Either `trimmomatic` or `fastp` can be used for read-quality trimming. Trimmomatic is used by default. Both tools trim low-quality regions of reads with a sliding window (with a window size of `trim_window_size`), cutting once the average quality within the window falls below `trim_quality_trim_score`. They will both discard the read if it is trimmed below `trim_min_length`. 

    By default, the trim_min_length is set to 75 bp. This is likely _too high_ for data generated using the Illumina VSP panel. We recommend setting this parameter to `50` in this case.

    If fastp is selected for analysis, fastp also implements the additional read-trimming parameters indicated below:

    | **Parameter** | **Explanation** |
    | --- | --- |
    | -g | enables polyG tail trimming |
    | -5 20 | enables read end-trimming |
    | -3 20 | enables read end-trimming |
    | --detect_adapter_for_pe | enables adapter-trimming **only for paired-end reads** |

    **Adapter removal**

    The `BBDuk` task removes adapters from sequence reads. To do this:

    - [Repair](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/repair-guide/) from the [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) package reorders reads in paired fastq files to ensure the forward and reverse reads of a pair are in the same position in the two fastq files.
    - [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/)  (*"Bestus Bioinformaticus" Decontamination Using Kmers*) is then used to trim the adapters and filter out all reads that have a 31-mer match to [PhiX](https://emea.illumina.com/products/by-type/sequencing-kits/cluster-gen-sequencing-reagents/phix-control-v3.html), which is commonly added to Illumina sequencing runs to monitor and/or improve overall run quality.
    
    ??? toggle "What are adapters and why do they need to be removed?"
        Adapters are manufactured oligonucleotide sequences attached to DNA fragments during the library preparation process. In Illumina sequencing, these adapter sequences are required for attaching reads to flow cells. You can read more about Illumina adapters [here](https://emea.support.illumina.com/bulletins/2020/06/illumina-adapter-portfolio.html). For genome analysis, it's important to remove these sequences since they're not actually from your sample. If you don't remove them, the downstream analysis may be affected.
        
    **Read Quantification**

    There are two methods for read quantification to choose from: [`fastq-scan`](https://github.com/rpetit3/fastq-scan) (default) or [`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Both quantify the forward and reverse reads in FASTQ files. In paired-end workflows, they also provide the total number of read pairs. This task is run once with raw reads as input and once with clean reads as input. If QC has been performed correctly, you should expect **fewer** clean reads than raw reads. `fastqc` also provides a graphical visualization of the read quality in an HTML file.

    **Read Identification (optional)**

    The `MIDAS` task is for the identification of reads to detect contamination with non-target taxa. This task is optional and turned off by default. It can be used by setting the `call_midas` input variable to `true`.

    The MIDAS reference database, located at **`gs://theiagen-large-public-files-rp/terra/theiaprok-files/midas/midas_db_v1.2.tar.gz`**, is provided as the default. It is possible to provide a custom database. More information is available [here](https://github.com/snayfach/MIDAS/blob/master/docs/ref_db.md).

    ??? toggle "How are the MIDAS output columns determined?"
        
        Example MIDAS report in the ****`midas_report` column:
        
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
  
    !!! techdetails "read_QC_trim Technical Details"
                
        |  | Links |
        | --- | --- |
        | Sub-workflow | [wf_read_QC_trim.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_read_QC_trim.wdl) |
        | Tasks | [task_fastp.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_fastp.wdl)<br>[task_trimmomatic.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_trimmomatic.wdl)<br>[task_bbduk.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_bbduk.wdl)<br>[task_fastq_scan.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_fastq_scan.wdl)<br>[task_midas.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/contamination/task_midas.wdl) |
        | Software Source Code | [fastp](https://github.com/OpenGene/fastp); [Trimmomatic](https://github.com/usadellab/Trimmomatic); [fastq-scan](https://github.com/rpetit3/fastq-scan); [MIDAS](https://github.com/snayfach/MIDAS)|
        | Software Documentation | [fastp](https://github.com/OpenGene/fastp); [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic); [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/); [fastq-scan](https://github.com/rpetit3/fastq-scan); [MIDAS](https://github.com/snayfach/MIDAS) |
        | Original Publication(s) | [Trimmomatic: a flexible trimmer for Illumina sequence data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4103590/)<br>[fastp: an ultra-fast all-in-one FASTQ preprocessor](https://academic.oup.com/bioinformatics/article/34/17/i884/5093234?login=false)<br>[An integrated metagenomics pipeline for strain profiling reveals novel patterns of bacterial transmission and biogeography](https://pubmed.ncbi.nlm.nih.gov/27803195/) |

??? task "`kraken2`: Taxonomic Classification"

    Kraken2 is a bioinformatics tool originally designed for metagenomic applications. It has additionally proven valuable for validating taxonomic assignments and checking contamination of single-species (e.g. bacterial isolate, eukaryotic isolate, viral isolate, etc.) whole genome sequence data.

    Kraken2 is run on the clean reads that result from the `read_QC_trim` subworkflow. By default, the Kraken2 database is set to the `k2_viral_20240112` database, located at `"gs://theiagen-large-public-files-rp/terra/databases/kraken2/k2_viral_20240112.tar.gz"`.

    !!! info "Database-dependent"
        The Kraken2 software is database-dependent and **taxonomic assignments are highly sensitive to the database used**. An appropriate database should contain the expected organism(s) (e.g. _Escherichia coli_) and other taxa that may be present in the reads (e.g. _Citrobacter freundii_, a common contaminant).

    !!! techdetails "Kraken2 Technical Details"    
        |  | Links |
        | --- | --- |
        | Task | [task_kraken2.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/contamination/task_kraken2.wdl) |
        | Software Source Code | [Kraken2 on GitHub](https://github.com/DerrickWood/kraken2/) |
        | Software Documentation | <https://github.com/DerrickWood/kraken2/wiki> |
        | Original Publication(s) | [Improved metagenomic analysis with Kraken 2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0) |

??? task "`KrakenTools extract_kraken_reads`: Read Binning"
    KrakenTools is a collection of scripts that can be used to help downstream analysis of Kraken2 results. In particular, this task uses the `extract_kraken_reads` script, which extracts reads classified at any user-specified taxonomy IDs. All parent and children reads of the specified taxonomic ID are also extracted.

    !!! techdetails "KrakenTools Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_kraken_tools.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/task_krakentools.wdl)
        | Software Source Code | [KrakenTools on GitHub](https://github.com/jenniferlu717/KrakenTools) |
        | Software Documentation | [KrakenTools on GitHub](https://github.com/jenniferlu717/KrakenTools) |
        | Original Publication | [Metagenome analysis using the Kraken software suite](https://doi.org/10.1038/s41596-022-00738-y) |

### Outputs

<div class="searchable-table" markdown="1">

| **Variable** | **Type** | **Description** |
|---|---|---|
| identified_organisms | Array[String] | A list of organisms that were able to be identified in the sample with the specified Kraken2 database |
| kraken2_classified_report | File | Standard Kraken2 output report. TXT filetype, but can be opened in Excel as a TSV file |
| kraken2_database | String | The name of the database used to run Kraken2 |
| kraken2_docker | String | Docker image used to run kraken2 |
| kraken2_report | File | Text document describing taxonomic prediction of every FASTQ record. This file can be very large and cumbersome to open and view |
| kraken2_version | String | The version of Kraken2 used in the analysis |
| results_by_taxon_tsv | File | A TSV file that contains the results for every taxon ID provided in the taxon_ids input variable that had reads identified; characterization (if applicable) and basic statistics regarding read count, assembly generation (if applicable), and general quality, are also associated with each bin |
| theiameta_panel_illumina_pe_analysis_date | String | Date the workflow was run |
| theiameta_panel_illumina_pe_version | String | Version of PHB used to run the workflow |

</div>

## References (if applicable)

>**Trimmomatic:** Anthony M. Bolger and others, Trimmomatic: a flexible trimmer for Illumina sequence data, *Bioinformatics*, Volume 30, Issue 15, August 2014, Pages 2114–2120, <https://doi.org/10.1093/bioinformatics/btu170>
<!-- -->
>**Fastq-Scan:** <https://github.com/rpetit3/fastq-scan>
<!-- -->
>**metaSPAdes:** Sergey Nurk and others, metaSPAdes: a new versatile metagenomic assembler, *Genome Res.* 2017 May; 27(5): 824–834., <https://doi.org/10.1101%2Fgr.213959.116>
<!-- -->
>**Pilon:** Bruce J. Walker and others. Pilon: An Integrated Tool for Comprehensive Microbial Variant Detection and Genome Assembly Improvement. *Plos One.* November 19, 2014. <https://doi.org/10.1371/journal.pone.0112963>
<!-- -->
>**Minimap2:** Heng Li, Minimap2: pairwise alignment for nucleotide sequences, *Bioinformatics*, Volume 34, Issue 18, September 2018, Pages 3094–3100, <https://doi.org/10.1093/bioinformatics/bty191>
<!-- -->
>**QUAST:** Alexey Gurevich and others, QUAST: quality assessment tool for genome assemblies, *Bioinformatics*, Volume 29, Issue 8, April 2013, Pages 1072–1075, <https://doi.org/10.1093/bioinformatics/btt086>
<!-- -->
>**Samtools:** Li, Heng, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils Homer, Gabor Marth, Goncalo Abecasis, Richard Durbin, and 1000 Genome Project Data Processing Subgroup. 2009. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25(16): 2078-2079.
<!-- -->
>**Bcftools:** Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li. Twelve years of SAMtools and BCFtools. GigaScience, Volume 10, Issue 2, February 2021, giab008, <https://doi.org/10.1093/gigascience/giab008>
<!-- -->
