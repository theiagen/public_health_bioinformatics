# TheiaEuk 

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibliity** | **Workflow Level** |
|---|---|---|---|---|
| [Genomic Characterization](../../workflows_overview/workflows_type.md/#genomic-characterization) | [Mycotics](../../workflows_overview/workflows_kingdom.md/#mycotics) | PHB v2.2.0 | Yes | Sample-level |

## TheiaEuk Workflows

**The TheiaEuk_PE workflow is for the assembly, quality assessment, and characterization of fungal genomes.** It is designed to accept Illumina paired-end sequencing data as the primary input. **It is currently intended only for haploid fungal genomes like _Candida auris_.** Analyzing diploid genomes using TheiaEuk should be attempted only with expert attention to the resulting genome quality.

All input reads are processed through "core tasks" in each workflow. The core tasks include raw-read quality assessment, read cleaning (quality trimming and adapter removal), de novo assembly, assembly quality assessment, and species taxon identification. For some taxa identified, "taxa-specific sub-workflows" will be automatically activated, undertaking additional taxa-specific characterization steps, including clade-typing and/or antifungal resistance detection.

!!! caption "TheiaEuk Workflow Diagram"
    ![TheiaEuk Workflow Diagram](../../assets/figures/TheiaEuk_Illumina_PE.png){width=75%}

### Inputs

!!! info "Input read data"

    The TheiaEuk_PE workflow takes in Illumina paired-end read data. Read file names should end with `.fastq` or `.fq`, with the optional addition of `.gz`. When possible, Theiagen recommends zipping files with [gzip](https://www.gnu.org/software/gzip/) prior to Terra upload to minimize data upload time.

    By default, the workflow anticipates 2 x 150bp reads (i.e. the input reads were generated using a 300-cycle sequencing kit). Modifications to the optional parameter for `trim_minlen` may be required to accommodate shorter read data, such as the 2 x 75bp reads generated using a 150-cycle sequencing kit.

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| theiaeuk_pe | **read1** | File | Unprocessed Illumina forward read file |  | Required |
| theiaeuk_pe | **read2** | File | Unprocessed Illumina reverse read file |  | Required |
| theiaeuk_pe | **samplename** | String | Name of Terra datatable |  | Required |
| busco | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| busco | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| busco | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/ezlabgva/busco:v5.3.2_cv1 | Optional |
| cg_pipeline_clean | **cg_pipe_opts** | String | Options to pass to CG-Pipeline for clean read assessment  | --fast | Optional |
| cg_pipeline_clean | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| cg_pipeline_clean | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/lyveset:1.1.4f | Optional |
| cg_pipeline_raw | **cg_pipe_opts** | String | Options to pass to CG-Pipeline for clean read assessment  | --fast | Optional |
| cg_pipeline_raw | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| cg_pipeline_raw | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/lyveset:1.1.4f | Optional |
| clean_check_reads | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| clean_check_reads | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| clean_check_reads | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/bactopia/gather_samples:2.0.2 | Optional |
| clean_check_reads | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional |
| clean_check_reads | **organism** | String | Internal component, do not modify |  | Do Not Modify, Optional |
| clean_check_reads | **workflow_series** | String | Internal component, do not modify |  | Do Not Modify, Optional |
| gambit | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| gambit | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/gambit:0.5.0 | Optional |
| merlin_magic | **agrvate_docker_image** | String | Internal component, do not modify | "us-docker.pkg.dev/general-theiagen/biocontainers/agrvate:1.0.2--hdfd78af_0" | Do Not Modify, Optional |
| merlin_magic | **assembly_only** | Boolean | Internal component, do not modify |  | Do Not Modify, Optional |
| merlin_magic | **call_poppunk** | Boolean | Internal component, do not modify | TRUE | Do Not Modify, Optional |
| merlin_magic | **call_shigeifinder_reads_input** | Boolean | Internal component, do not modify | FALSE | Do Not Modify, Optional |
| merlin_magic | **emmtypingtool_docker_image** | String | Internal component, do not modify | us-docker.pkg.dev/general-theiagen/staphb/emmtypingtool:0.0.1 | Do Not Modify, Optional |
| merlin_magic | **hicap_docker_image** | String | Internal component, do not modify | us-docker.pkg.dev/general-theiagen/biocontainers/hicap:1.0.3--py_0 | Do Not Modify, Optional |
| merlin_magic | **ont_data** | Boolean | Internal component, do not modify |  | Do Not Modify, Optional |
| merlin_magic | **paired_end** | Boolean | Internal component, do not modify |  | Do Not Modify, Optional |
| merlin_magic | **pasty_docker_image** | String | Internal component, do not modify | us-docker.pkg.dev/general-theiagen/staphb/pasty:1.0.3 | Do Not Modify, Optional |
| merlin_magic | **pasty_min_coverage** | Int | Internal component, do not modify | 95 | Do Not Modify, Optional |
| merlin_magic | **pasty_min_pident** | Int | Internal component, do not modify | 95 | Do Not Modify, Optional |
| merlin_magic | **shigatyper_docker_image** | String | Internal component, do not modify | us-docker.pkg.dev/general-theiagen/staphb/shigatyper:2.0.5 | Do Not Modify, Optional |
| merlin_magic | **shigeifinder_docker_image** | String | Internal component, do not modify | us-docker.pkg.dev/general-theiagen/staphb/shigeifinder:1.3.5 | Do Not Modify, Optional |
| merlin_magic | **snippy_query_gene** | String | Internal component, do not modify |  | Do Not Modify, Optional |
| merlin_magic | **srst2_gene_max_mismatch** | Int | Internal component, do not modify | 2000 | Do Not Modify, Optional |
| merlin_magic | **srst2_max_divergence** | Int | Internal component, do not modify | 20 | Do Not Modify, Optional |
| merlin_magic | **srst2_min_cov** | Int | Internal component, do not modify | 80 | Do Not Modify, Optional |
| merlin_magic | **srst2_min_depth** | Int | Internal component, do not modify | 5 | Do Not Modify, Optional |
| merlin_magic | **srst2_min_edge_depth** | Int | Internal component, do not modify | 2 | Do Not Modify, Optional |
| merlin_magic | **staphopia_sccmec_docker_image** | String | Internal component, do not modify | us-docker.pkg.dev/general-theiagen/biocontainers/staphopia-sccmec:1.0.0--hdfd78af_0 | Do Not Modify, Optional |
| merlin_magic | **tbp_parser_coverage_threshold** | Int | Internal component, do not modify | 100 | Do Not Modify, Optional |
| merlin_magic | **tbp_parser_debug** | Boolean | Internal component, do not modify | FALSE | Do Not Modify, Optional |
| merlin_magic | **tbp_parser_docker_image** | String | Internal component, do not modify | us-docker.pkg.dev/general-theiagen/theiagen/tbp-parser:1.3.6 | Do Not Modify, Optional |
| merlin_magic | **tbp_parser_min_depth** | Int | Internal component, do not modify | 10 | Do Not Modify, Optional |
| merlin_magic | **tbp_parser_operator** | String | Internal component, do not modify | "Operator not provided" | Do Not Modify, Optional |
| merlin_magic | **tbp_parser_output_seq_method_type** | String | Internal component, do not modify | "WGS" | Do Not Modify, Optional |
| merlin_magic | **tbp_parser_output_seq_method_type** | String | Internal component, do not modify | "Sequencing method not provided" | Do Not Modify, Optional |
| merlin_magic | **tbprofiler_additional_outputs** | Boolean | Internal component, do not modify | FALSE | Do Not Modify, Optional |
| merlin_magic | **tbprofiler_cov_frac_threshold** | Int | Internal component, do not modify | 1 | Do Not Modify, Optional |
| merlin_magic | **tbprofiler_custom_db** | File | Internal component, do not modify |  | Do Not Modify, Optional |
| merlin_magic | **tbprofiler_mapper** | String | Internal component, do not modify | bwa | Do Not Modify, Optional |
| merlin_magic | **tbprofiler_min_af** | Float | Internal component, do not modify | 0.1 | Do Not Modify, Optional |
| merlin_magic | **tbprofiler_min_af_pred** | Float | Internal component, do not modify | 0.1 | Do Not Modify, Optional |
| merlin_magic | **tbprofiler_min_depth** | Int | Internal component, do not modify | 10 | Do Not Modify, Optional |
| merlin_magic | **tbprofiler_run_custom_db** | Boolean | Internal component, do not modify | FALSE | Do Not Modify, Optional |
| merlin_magic | **tbprofiler_variant_caller** | String | Internal component, do not modify | freebayes | Do Not Modify, Optional |
| merlin_magic | **tbprofiler_variant_calling_params** | String | Internal component, do not modify | None | Do Not Modify, Optional |
| merlin_magic | **virulencefinder_coverage_threshold** | Float | Internal component, do not modify |  | Do Not Modify, Optional |
| merlin_magic | **virulencefinder_database** | String | Internal component, do not modify | "virulence_ecoli" | Do Not Modify, Optional |
| merlin_magic | **virulencefinder_docker_image** | String | Internal component, do not modify | us-docker.pkg.dev/general-theiagen/staphb/virulencefinder:2.0.4 | Do Not Modify, Optional |
| merlin_magic | **virulencefinder_identity_threshold** | Float | Internal component, do not modify |  | Do Not Modify, Optional |
| qc_check_task | **ani_highest_percent** | Float | Internal component, do not modify |  | Do Not Modify, Optional |
| qc_check_task | **ani_highest_percent_bases_aligned** | Float | Internal component, do not modify |  | Do Not Modify, Optional |
| qc_check_task | **assembly_length_unambiguous** | Int | Internal component, do not modify |  | Do Not Modify, Optional |
| qc_check_task | **assembly_mean_coverage** | Float | Internal component, do not modify |  | Do Not Modify, Optional |
| qc_check_task | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional |
| qc_check_task | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| qc_check_task | **docker** | String | The Docker container to use for the task |  "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16" | Optional |
| qc_check_task | **kraken_human** | Float | Internal component, do not modify |  | Do Not Modify, Optional |
| qc_check_task | **kraken_human_dehosted** | Float | Internal component, do not modify |  | Do Not Modify, Optional |
| qc_check_task | **kraken_sc2** | Float | Internal component, do not modify |  | Do Not Modify, Optional |
| qc_check_task | **kraken_sc2_dehosted** | Float | Internal component, do not modify |  | Do Not Modify, Optional |
| qc_check_task | **kraken_target_organism** | Float | Internal component, do not modify |  | Do Not Modify, Optional |
| qc_check_task | **kraken_target_organism_dehosted** | Float | Internal component, do not modify |  | Do Not Modify, Optional |
| qc_check_task | **meanbaseq_trim** | String | Internal component, do not modify |  | Do Not Modify, Optional |
| qc_check_task | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| qc_check_task | **midas_secondary_genus_abundance** | Int | Internal component, do not modify |  | Do Not Modify, Optional |
| qc_check_task | **midas_secondary_genus_coverage** | Float | Internal component, do not modify |  | Do Not Modify, Optional |
| qc_check_task | **number_Degenerate** | Int | Internal component, do not modify |  | Do Not Modify, Optional |
| qc_check_task | **number_N** | Int | Internal component, do not modify |  | Do Not Modify, Optional |
| qc_check_task | **percent_reference_coverage** | Float | Internal component, do not modify |  | Do Not Modify, Optional |
| qc_check_task | **sc2_s_gene_mean_coverage** | Float | Internal component, do not modify |  | Do Not Modify, Optional |
| qc_check_task | **sc2_s_gene_percent_coverage** | Float | Internal component, do not modify |  | Do Not Modify, Optional |
| qc_check_task | **vadr_num_alerts** | String | Internal component, do not modify |  | Do Not Modify, Optional |
| quast | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| quast | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/quast:5.0.2 | Optional |
| quast | **min_contig_length** | Int | Minimum length of contig for QUAST | 500 | Optional |
| rasusa_task | **bases** | String | Explicitly set the number of bases required e.g., 4.3kb, 7Tb, 9000, 4.1MB |  | Optional |
| rasusa_task | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional |
| rasusa_task | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| rasusa_task | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/rasusa:0.7.0 | Optional |
| rasusa_task | **frac** | Float | Subsample to a fraction of the reads |  | Optional |
| rasusa_task | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| rasusa_task | **num** | Int | Subsample to a specific number of reads |  | Optional |
| rasusa_task | **seed** | Int | Random seed to use |  | Optional |
| raw_check_reads | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| raw_check_reads | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| raw_check_reads | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/bactopia/gather_samples:2.0.2 | Optional |
| raw_check_reads | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional |
| raw_check_reads | **organism** | String | Internal component, do not modify |  | Do Not Modify, Optional |
| raw_check_reads | **workflow_series** | String | Internal component, do not modify |  | Do Not Modify, Optional |
| read_QC_trim | **adapters** | File | File with adapter sequences to be removed |  | Optional |
| read_QC_trim | **bbduk_mem** | Int | Memory allocated to the BBDuk VM | 8 | Optional |
| read_QC_trim | **call_kraken** | Boolean | If true, Kraken2 is executed on the dataset | FALSE | Optional |
| read_QC_trim | **call_midas** | Boolean | Internal component, do not modify | FALSE | Do Not Modify, Optional |
| read_QC_trim | **fastp_args** | String | Additional arguments to pass to fastp | --detect_adapter_for_pe -g -5 20 -3 20 | Optional |
| read_QC_trim | **kraken_db** | File | Database to use with kraken2 |  | Optional |
| read_QC_trim | **kraken_disk_size** | Int | Amount of storage (in GB) to allocate to the task |  | Optional |
| read_QC_trim | **kraken_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task |  | Optional |
| read_QC_trim | **midas_db** | File | Internal component, do not modify |  | Do Not Modify, Optional |
| read_QC_trim | **phix** | File | A file containing the phix used during Illumina sequencing; used in the BBDuk task |  | Optional |
| read_QC_trim | **read_processing** | String | Read trimming software to use, either "trimmomatic" or "fastp" | trimmomatic | Optional |
| read_QC_trim | **read_qc** | String | Allows the user to decide between fastq_scan (default) and fastqc for the evaluation of read quality. | "fastq_scan" | Optional |
| read_QC_trim | **target_organism** | String | This string is searched for in the kraken2 outputs to extract the read percentage |  | Optional |
| read_QC_trim | **trim_minlength** | Int | Specifies minimum length of each read after trimming to be kept | 75 | Optional |
| read_QC_trim | **trim_quality_trim_score** | Int | Specifies the average quality of bases in a sliding window to be kept | 20 | Optional |
| read_QC_trim | **trim_window_size** | Int | Specifies window size for trimming (the number of bases to average the quality across) | 10 | Optional |
| read_QC_trim | **trimmomatic_args** | String | Additional arguments for trimmomatic |  | Optional |
| read_QC_trim | **workflow_series** | String | Internal component, do not modify |  | Do Not Modify, Optional |
| shovill_pe | **assembler** | String | Assembler to use (spades, skesa, velvet or megahit), see <https://github.com/tseemann/shovill#--assembler> | "skesa" | Optional |
| shovill_pe | **assembler_options** | String | Assembler-specific options that you might choose, see <https://github.com/tseemann/shovill#--opts> |  | Optional |
| shovill_pe | **depth** | Int | User specified depth of coverage for downsampling (see <https://github.com/tseemann/shovill#--depth and https://github.com/tseemann/shovill#main-steps>) | 150 | Optional |
| shovill_pe | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| shovill_pe | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/shovill:1.1.0 | Optional |
| shovill_pe | **genome_length** | String | Internal component, do not modify |  | Do Not Modify, Optional |
| shovill_pe | **kmers** | String | User-specified Kmer length to override choice made by Shovill, see <https://github.com/tseemann/shovill#--kmers> | auto | Optional |
| shovill_pe | **min_contig_length** | Int | Minimum contig length to keep in final assembly  | 200 | Optional |
| shovill_pe | **min_coverage** | Float | Minimum contig coverage to keep in final assembly | 2 | Optional |
| shovill_pe | **nocorr** | Boolean | Disable correction of minor assembly errors by Shovill (see <https://github.com/tseemann/shovill#main-steps>) | FALSE | Optional |
| shovill_pe | **noreadcorr** | Boolean | Disable correction of sequencing errors in reads by Shovill (see <https://github.com/tseemann/shovill#main-steps>) | FALSE | Optional |
| shovill_pe | **nostitch** | Boolean | Disable read stitching by Shovill (see <https://github.com/tseemann/shovill#main-steps>) | FALSE | Optional |
| shovill_pe | **trim** | Boolean | Enable adaptor trimming (see <https://github.com/tseemann/shovill#main-step>s) | FALSE | Optional |
| theiaeuk_pe | **busco_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| theiaeuk_pe | **call_rasusa** | Boolean | If true, launch rasusa task to subsample raw reads to read depth of 150X | TRUE | Optional |
| theiaeuk_pe | **gambit_db_genomes** | File | User-provided database of assembled query genomes; requires complementary signatures file. If not provided, uses default database, "/gambit-db" | gs://gambit-databases-rp/1.3.0/gambit-metadata-1.3-231016.gdb | Optional |
| theiaeuk_pe | **gambit_db_signatures** | File | User-provided signatures file; requires complementary genomes file. If not specified, the file from the docker container will be used.  | gs://gambit-databases-rp/1.3.0/gambit-signatures-1.3-231016.gs | Optional |
| theiaeuk_pe | **genome_length** | Int | User-specified expected genome size to be used in genome statistics calculations |  | Optional |
| theiaeuk_pe | **max_genome_size** | Int | Maximum genome size able to pass read screening | 50000000 | Optional |
| theiaeuk_pe | **min_basepairs** | Int | Minimum number of base pairs able to pass read screening | 2241820 | Optional |
| theiaeuk_pe | **min_coverage** | Int | Minimum genome coverage able to pass read screening | 10 | Optional |
| theiaeuk_pe | **min_genome_size** | Int | Minimum genome size able to pass read screening | 100000 | Optional |
| theiaeuk_pe | **min_proportion** | Int | Minimum proportion of total reads in each read file to pass read screening | 50 | Optional |
| theiaeuk_pe | **min_reads** | Int | Minimum number of reads to pass read screening | 10000 | Optional |
| theiaeuk_pe | **skip_screen** | Boolean | Option to skip the read screening prior to analysis | FALSE | Optional |
| theiaeuk_pe | **subsample_coverage** | Float | Read depth for RASUSA task to subsample reads to | 150 | Optional |
| version_capture | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional |

### Workflow tasks (performed for all taxa)

??? task "`versioning`: Version capture for TheiaEuk"

    The `versioning` task captures the workflow version from the GitHub (code repository) version.
        
    !!! techdetails "Version Capture Technical details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_versioning.wdl](https://github.com/theiagen/public_health_bacterial_genomics/blob/main/tasks/task_versioning.wdl) |

??? task "`screen`: Total Raw Read Quantification and Genome Size Estimation"

    The [`screen`](https://github.com/theiagen/public_health_bacterial_genomics/blob/main/tasks/quality_control/task_screen.wdl) task ensures the quantity of sequence data is sufficient to undertake genomic analysis. It uses bash commands for quantification of reads and base pairs, and [mash](https://mash.readthedocs.io/en/latest/index.html) sketching to estimate the genome size and its coverage. At each step, the results are assessed relative to pass/fail criteria and thresholds that may be defined by optional user inputs. Samples that do not meet these criteria will not be processed further by the workflow:

    1. Total number of reads: A sample will fail the read screening task if its total number of reads is less than or equal to `min_reads`.
    2. The proportion of basepairs reads in the forward and reverse read files: A sample will fail the read screening if fewer than `min_proportion` basepairs are in either the reads1 or read2 files.
    3. Number of basepairs: A sample will fail the read screening if there are fewer than `min_basepairs` basepairs
    4. Estimated genome size:  A sample will fail the read screening if the estimated genome size is smaller than `min_genome_size` or bigger than `max_genome_size`.
    5. Estimated genome coverage: A sample will fail the read screening if the estimated genome coverage is less than the `min_coverage`.

    Read screening is undertaken on both the raw and cleaned reads. The task may be skipped by setting the `skip_screen` variable to true.

    Default values vary between the PE and SE workflow. The rationale for these default values can be found below.
    
    | Variable  | Rationale |
    | --- | --- |
    | `skip_screen` | Prevent the read screen from running |
    | `min_reads` | Minimum number of base pairs for 20x coverage of _Hansenula polymorpha_  divided by 300 (longest Illumina read length) |
    | `min_basepairs` | Greater than 10x coverage of _Hansenula polymorpha_  |
    | `min_genome_size` | Based on the _Hansenula polymorpha_  genome - the smallest fungal genome as of 2015-04-02 (8.97 Mbp) |
    | `max_genome_size` | Based on the _Cenococcum geophilum_  genome, the biggest pathogenic fungal genome, (177.57 Mbp) |
    | `min_coverage` | A bare-minimum coverage for genome characterization. Higher coverage would be required for high-quality phylogenetics. |
    | `min_proportion` | Greater than 50% reads are in the read1 file; others are in the read2 file |

    !!! techdetails "Screen Technical Details"
        
        There is a single WDL task for read screening. The `screen` task is run twice, once for raw reads and once for clean reads.
        
        |  | Links |
        | --- | --- |
        | Task | [task_screen.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/task_screen.wdl)  |

??? task "`rasusa`: Read subsampling"

    The RASUSA task performs subsampling of the raw reads. By default, this task will subsample reads to a depth of 150X using the estimated genome length produced during the preceding raw read screen. The user can prevent the task from being launched by setting the `call_rasusa`variable to false. 

    The user can also provide an estimated genome length for the task to use for subsampling using the `genome_size` variable. In addition, the read depth can be modified using the `subsample_coverage` variable.
        
    !!! techdetails "RASUSA Technical Details"

        |  | TheiaEuk_Illumina_PE_PHB |
        | --- | --- |
        | Task | [task_rasusa.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/task_rasusa.wdl) |

??? task "`read_QC_trim`: Read Quality Trimming, Adapter Removal, Quantification, and Identification"

    `read_QC_trim` is a sub-workflow within TheiaEuk that removes low-quality reads, low-quality regions of reads, and sequencing adapters to improve data quality. It uses a number of tasks, described below.

    **Read quality trimming**

    Either `trimmomatic` or `fastp` can be used for read-quality trimming. Trimmomatic is used by default. Both tools trim low-quality regions of reads with a sliding window (with a window size of `trim_window_size`), cutting once the average quality within the window falls below `trim_quality_trim_score`. They will both discard the read if it is trimmed below `trim_minlen`. 

    If fastp is selected for analysis, fastp also implements the additional read-trimming steps indicated below:

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

    There are two methods for read quantification to choose from: [`fastq-scan`](https://github.com/rpetit3/fastq-scan) (default) or [`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Both quantify the forward and reverse reads in FASTQ files. In TheiaProk_Illumina_PE, they also provide the total number of read pairs. This task is run once with raw reads as input and once with clean reads as input. If QC has been performed correctly, you should expect **fewer** clean reads than raw reads. `fastqc` also provides a graphical visualization of the read quality.

    **Read Identification (optional)**

    The `MIDAS` task is for the identification of reads to detect contamination with non-target taxa. This task is optional and turned off by default. It can be used by setting the `call_midas` input variable to `true`.

    The MIDAS tool was originally designed for metagenomic sequencing data but has been co-opted for use with bacterial isolate WGS methods. It can be used to detect contamination present in raw sequencing data by estimating bacterial species abundance in bacterial isolate WGS data. If a secondary genus is detected above a relative frequency of 0.01 (1%), then the sample should fail QC and be investigated further for potential contamination.

    This task is similar to those used in commercial software, BioNumerics, for estimating secondary species abundance.

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
        
        The value in the `midas_primary_genus` column is derived by ordering the rows in order of "relative_abundance" and identifying the genus of top species in the "species_id" column (Salmonella). The value in the `midas_secondary_genus` column is derived from the genus of the second-most prevalent genus in the "species_id" column (Citrobacter). The `midas_secondary_genus_abundance` column is the "relative_abundance" of the second-most prevalent genus (0.009477003). The `midas_secondary_genus_coverage` is the "coverage" of the second-most prevalent genus (0.995216227).
        
    !!! techdetails "read_QC_trim Technical Details"
                
        |  | Links |
        | --- | --- |
        | Sub-workflow | [wf_read_QC_trim.wdl](https://github.com/theiagen/public_health_bacterial_genomics/blob/main/workflows/wf_read_QC_trim.wdl) |
        | Tasks | [task_fastp.wdl](https://github.com/theiagen/public_health_bacterial_genomics/blob/main/tasks/quality_control/task_fastp.wdl)<br>[task_trimmomatic.wdl](https://github.com/theiagen/public_health_bacterial_genomics/blob/main/tasks/quality_control/task_trimmomatic.wdl#L3) (PE subtask)<br>[task_bbduk.wdl](https://github.com/theiagen/public_health_bacterial_genomics/blob/main/tasks/quality_control/task_bbduk.wdl)<br>[task_fastq_scan.wdl](https://github.com/theiagen/public_health_bacterial_genomics/blob/main/tasks/quality_control/task_fastq_scan.wdl#L3) (PE subtask)<br>[task_midas.wdl](https://github.com/theiagen/public_health_bacterial_genomics/blob/main/tasks/taxon_id/task_midas.wdl)<br>[task_kraken2.wdl](https://github.com/theiagen/public_health_bacterial_genomics/blob/main/tasks/taxon_id/task_kraken2.wdl) |
        | Software Source Code | [fastp](https://github.com/OpenGene/fastp); [Trimmomatic](https://github.com/usadellab/Trimmomatic); [fastq-scan](https://github.com/rpetit3/fastq-scan); [MIDAS](https://github.com/snayfach/MIDAS); [Kraken2](https://github.com/DerrickWood/kraken2)|
        | Software Documentation | [fastp](https://github.com/OpenGene/fastp); [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic); [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/); [fastq-scan](https://github.com/rpetit3/fastq-scan); [MIDAS](https://github.com/snayfach/MIDAS); [Kraken2](https://github.com/DerrickWood/kraken2/wiki) |
        | Original Publication(s) | *[Trimmomatic: a flexible trimmer for Illumina sequence data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4103590/)<br>*[fastp: an ultra-fast all-in-one FASTQ preprocessor](https://academic.oup.com/bioinformatics/article/34/17/i884/5093234?login=false)<br>*[An integrated metagenomics pipeline for strain profiling reveals novel patterns of bacterial transmission and biogeography](https://pubmed.ncbi.nlm.nih.gov/27803195/)<br>*[Improved metagenomic analysis with Kraken 2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0) |

??? task "`shovill`: _De novo_ Assembly"

    De Novo assembly will be undertaken only for samples that have sufficient read quantity and quality, as determined by the `screen` task assessment of clean reads. 

    In TheiaEuk, assembly is performed using the [Shovill](https://github.com/tseemann/shovill) pipeline. This undertakes the assembly with one of four assemblers ([SKESA](https://github.com/ncbi/SKESA) (default), [SPAdes](https://github.com/ablab/spades), [Velvet](https://github.com/dzerbino/velvet/), [Megahit](https://github.com/voutcn/megahit)), but also performs [a number of pre- and post-processing steps](https://github.com/tseemann/shovill#main-steps) to improve the resulting genome assembly. Shovill uses an estimated genome size (see [here](https://github.com/tseemann/shovill#--gsize)). If this is not provided by the user as an optional input, Shovill will estimate the genome size using [mash](https://mash.readthedocs.io/en/latest/index.html). Adaptor trimming can be undertaken with Shovill by setting the `trim` option to "true", but this is set to "false" by default as [alternative adapter trimming](https://www.notion.so/TheiaProk-Workflow-Series-89b9c08406094ec78d08a578fe861626?pvs=21) is undertaken in the TheiaEuk workflow.

    ??? toggle "What is _de novo_  assembly?"
        _De novo_  assembly is the process or product of attempting to reconstruct a genome from scratch (without prior knowledge of the genome) using sequence reads. Assembly of fungal genomes from short-reads will produce multiple contigs per chromosome rather than a single contiguous sequence for each chromosome.
        
    !!! techdetails "Shovill Technical Details"
        |  | Links |
        | --- | --- |
        | TheiaEuk WDL Task | [task_shovill.wdl](https://github.com/theiagen/public_health_bacterial_genomics/blob/main/tasks/assembly/task_shovill.wdl#L3) |
        | Software code repository and documentation | [Shovill on GitHub](https://github.com/tseemann/shovill) |

??? task "`QUAST`: Assembly Quality Assessment"

    [`QUAST`](https://github.com/ablab/quast) (**QU**ality **AS**sessment **T**ool) evaluates genome assemblies by computing several metrics that describe the assembly quality, including the total number of bases in the assembly, the length of the largest contig in the assembly, and the assembly percentage GC content.

    !!! techdetails "QUAST Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_quast.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/task_quast.wdl) |
        | Software Source Code | [QUAST on GitHub](https://github.com/ablab/quast) |
        | Software Documentation | https://quast.sourceforge.net/docs/manual.html |
        | Orginal publication | [QUAST: quality assessment tool for genome assemblies](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3624806/) |

??? task "`CG-Pipeline`: Assessment of Read Quality, and Estimation of Genome Coverage"

    The`cg_pipeline` task generates metrics about read quality and estimates the coverage of the genome using the "run_assembly_readMetrics.pl" script from [CG-Pipeline](https://github.com/lskatz/CG-Pipeline/). The genome coverage estimates are calculated using both using raw and cleaned reads, using either a user-provided `genome_size` or the estimated genome length generated by QUAST.

    !!! techdetails "CG-Pipeline Technical Details"
        The `cg_pipeline` task is run twice in TheiaEuk, once with raw reads, and once with clean reads.
        
        |  | Links |
        | --- | --- |
        | Task | [task_cg_pipeline.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/task_cg_pipeline.wdl) |
        | Software Source Code | [CG-Pipeline on GitHub](https://github.com/lskatz/CG-Pipeline/) |
        | Software Documentation | [CG-Pipeline on GitHub](https://github.com/lskatz/CG-Pipeline/) |
        | Original Publication(s) | [A computational genomics pipeline for prokaryotic sequencing projects](https://academic.oup.com/bioinformatics/article/26/15/1819/188418) |

??? task "`GAMBIT`: **Taxon Assignment**"

    [`GAMBIT`](https://github.com/jlumpe/gambit) determines the taxon of the genome assembly using a k-mer based approach to match the assembly sequence to the closest complete genome in a database, thereby predicting its identity. Sometimes, GAMBIT can confidently designate the organism to the species level. Other times, it is more conservative and assigns it to a higher taxonomic rank.

    For additional details regarding the GAMBIT tool and a list of available GAMBIT databases for analysis, please consult the [GAMBIT](https://www.notion.so/GAMBIT-7c1376b861d0486abfbc316480046bdc?pvs=21) tool documentation.

    !!! techdetails "GAMBIT Technical Details"

        |  | Links |
        | --- | --- |
        | Task | [task_gambit.wdl](https://github.com/theiagen/public_health_bacterial_genomics/blob/main/tasks/taxon_id/task_gambit.wdl) |
        | Software Source Code | [GAMBIT on GitHub](https://github.com/jlumpe/gambit) |
        | Software Documentation | [GAMBIT ReadTheDocs](https://gambit-genomics.readthedocs.io/en/latest/) |
        | Original Publication(s) | [GAMBIT (Genomic Approximation Method for Bacterial Identification and Tracking): A methodology to rapidly leverage whole genome sequencing of bacterial isolates for clinical identification](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0277575) |

??? task "**`QC_check`: Check QC Metrics Against User-Defined Thresholds (optional)**"

    The `qc_check` task compares generated QC metrics against user-defined thresholds for each metric. This task will run if the user provides a `qc_check_table` .tsv file. If all QC metrics meet the threshold, the `qc_check` output variable will read `QC_PASS`. Otherwise, the output will read `QC_NA` if the task could not proceed or `QC_ALERT` followed by a string indicating what metric failed.

    The `qc_check` task applies quality thresholds according to the sample taxa. The sample taxa is taken from the `gambit_predicted_taxon` value inferred by the GAMBIT module OR can be manually provided by the user using the `expected_taxon` workflow input.

    ??? toggle "Formatting the _qc_check_table.tsv_"

        - The first column of the qc_check_table lists the taxa that the task will assess and the header of this column must be "taxon".
        - Any genus or species can be included as a row of the qc_check_table. However, these taxa must **uniquely** match the sample taxa, meaning that the file can include multiple species from the same genus (Vibrio_cholerae and Vibrio_vulnificus), but not both a genus row and species within that genus (Vibrio and Vibrio cholerae). **The taxa should be formatted with the first letter capitalized and underscores in lieu of spaces.**
        - Each subsequent column indicates a QC metric and lists a threshold for each taxa that will be checked. **The column names must exactly match expected values, so we highly recommend copy and pasting from the template files below.**

    ??? toggle "Template _qc_check_table.tsv_ files"

        TheiaEuk_Illumina_PE_PHB: [theiaeuk_qc_check_template.tsv](../../assets/files/TheiaEuk_qc_check_template.tsv)

        !!! warning "Example Purposes Only"
            QC threshold values shown are for example purposes only and should not be presumed to be sufficient for every dataset.
    
    !!! techdetails "QC_Check Technical Details"    
        
        |  | Links |
        | --- | --- |
        | Task | [task_qc_check.wdl](https://github.com/theiagen/public_health_bacterial_genomics/blob/main/tasks/quality_control/task_qc_check.wdl) |

### Organism-specific Characterization

The TheiaEuk workflow automatically activates taxa-specific tasks after identification of relevant taxa using `GAMBIT`. Many of these taxa-specific tasks do not require any additional workflow tasks from the user.

??? toggle "_Candida auris_"

    Two tools are deployed when _Candida auris is_  identified. First, the Cladetyping tool is launched to determine the clade of the specimen by comparing the sequence to five clade-specific reference files. The output of the clade typing task will be used to specify the reference genome for the antifungal resistance detection tool. To detect mutations that may confer antifungal resistance, `Snippy` is used to find all variants relative to the clade-specific reference, then these variants are queried for product names associated with resistance according to the MARDy database (<http://mardy.dide.ic.ac.uk/index.php>).

    **Default reference genomes used for clade typing and antimicrobial resistance gene detection of C. auris**

    | Clade | Genome Accession | Assembly Name | Strain | NCBI Submitter | Included mutations in AMR genes (not comprehensive) |
    | --- | --- | --- | --- | --- | --- |
    | Candida auris Clade I | GCA_002759435.2 | Cand_auris_B8441_V2 | B8441 | Centers for Disease Control and Prevention |  |
    | Candida auris Clade II | GCA_003013715.2 | ASM301371v2 | B11220 | Centers for Disease Control and Prevention |  |
    | Candida auris Clade III | GCA_002775015.1 | Cand_auris_B11221_V1 | B11221 | Centers for Disease Control and Prevention | _ERG11_ V125A/F126L |
    | Candida auris Clade IV | GCA_003014415.1 | Cand_auris_B11243 | B11243 | Centers for Disease Control and Prevention | _ERG11_ Y132F |
    | Candida auris Clade V | GCA_016809505.1 | ASM1680950v1 | IFRC2087 | Centers for Disease Control and Prevention |  |

    The genes in which there are known resistance-conferring mutations for this pathogen are:

    - FKS1
    - ERG11 (lanosterol 14-alpha demethylase)
    - FUR1 (uracil phosphoribosyltransferase)

    Mutations in these genes that are known to confer resistance are shown below (source: MARDy database http://mardy.dide.ic.ac.uk/index.php)

    | **Organism** | **Found in** | **Gene name** | **Gene locus** | **AA mutation** | **Drug** | **Tandem repeat name** | **Tandem repeat sequence** | **Reference** |
    | --- | --- | --- | --- | --- | --- | --- | --- | --- |
    | **Candida auris** | **Human** | **ERG11** |  | **Y132F** | **Fluconazole** |  |  | [**10.1093/cid/ciw691**](https://academic.oup.com/cid/article/64/2/134/2706620/Simultaneous-Emergence-of-Multidrug-Resistant) |
    | **Candida auris** | **Human** | **ERG11** |  | **K143R** | **Fluconazole** |  |  | [**10.1093/cid/ciw691**](https://academic.oup.com/cid/article/64/2/134/2706620/Simultaneous-Emergence-of-Multidrug-Resistant) |
    | **Candida auris** | **Human** | **ERG11** |  | **F126T** | **Fluconazole** |  |  | [**10.1093/cid/ciw691**](https://academic.oup.com/cid/article/64/2/134/2706620/Simultaneous-Emergence-of-Multidrug-Resistant) |
    | **Candida auris** | **Human** | **FKS1** |  | **S639P** | **Micafungin** |  |  | [**10.1016/j.diagmicrobio.2017.10.021**](https://www.sciencedirect.com/science/article/pii/S0732889317303498) |
    | **Candida auris** | **Human** | **FKS1** |  | **S639P** | **Caspofungin** |  |  | [**10.1016/j.diagmicrobio.2017.10.021**](https://www.sciencedirect.com/science/article/pii/S0732889317303498) |
    | **Candida auris** | **Human** | **FKS1** |  | **S639P** | **Anidulafungin** |  |  | [**10.1016/j.diagmicrobio.2017.10.021**](https://www.sciencedirect.com/science/article/pii/S0732889317303498) |
    | **Candida auris** | **Human** | **FKS1** |  | **S639F** | **Micafungin** |  |  | [**10.1093/jac/dkx480**](https://academic.oup.com/jac/advance-article/doi/10.1093/jac/dkx480/4794718) |
    | **Candida auris** | **Human** | **FKS1** |  | **S639F** | **Caspofungin** |  |  | [**10.1093/jac/dkx480**](https://academic.oup.com/jac/advance-article/doi/10.1093/jac/dkx480/4794718) |
    | **Candida auris** | **Human** | **FKS1** |  | **S639F** | **Anidulafungin** |  |  | [**10.1093/jac/dkx480**](https://academic.oup.com/jac/advance-article/doi/10.1093/jac/dkx480/4794718) |
    | **Candida auris** | **Human** | **FUR1** | **CAMJ_004922** | **F211I** | **5-flucytosine** |  |  | [**https://doi.org/10.1038/s41426-018-0045-x**](https://www.nature.com/articles/s41426-018-0045-x) |

??? toggle "_Candida albicans_"

    When this species is detected by the taxon ID tool, an antifungal resistance detection task is deployed. To detect mutations that may confer antifungal resistance, `Snippy` is used to find all variants relative to the clade-specific reference, and these variants are queried for product names associated with resistance according to the MARDy database (<http://mardy.dide.ic.ac.uk/index.php>).

    The genes in which there are known resistance-conferring mutations for this pathogen are:

    - ERG11
    - GCS1 (FKS1)
    - FUR1
    - RTA2

??? toggle "_Aspergillus fumigatus_"

    When this species is detected by the taxon ID tool an antifungal resistance detection task is deployed. To detect mutations that may confer antifungal resistance, `Snippy` is used to find all variants relative to the clade-specific reference, and these variants are queried for product names associated with resistance according to the MARDy database (<http://mardy.dide.ic.ac.uk/index.php>).

    The genes in which there are known resistance-conferring mutations for this pathogen are:

    - Cyp51A
    - HapE
    - COX10 (AFUA_4G08340)

??? toggle "_Cryptococcus neoformans_"

    When this species is detected by the taxon ID tool an antifungal resistance detection task is deployed. To detect mutations that may confer antifungal resistance, `Snippy` is used to find all variants relative to the clade-specific reference, and these variants are queried for product names associated with resistance according to the MARDy database (<http://mardy.dide.ic.ac.uk/index.php>).

    The gene in which there are known resistance-conferring mutations for this pathogen is:

    - ERG11 (CNA00300)

### Outputs

| **Variable** | **Type** | **Description** |
|---|---|---|
| cg_pipeline_docker | String | Docker file used for running CG-Pipeline on cleaned reads |
| cg_pipeline_report | File | TSV file of read metrics from raw reads, including average read length, number of reads, and estimated genome coverage |
| est_coverage_clean | Float | Estimated coverage calculated from   clean reads and genome length |
| est_coverage_raw | Float | Estimated coverage calculated from  raw reads and genome length |
| r1_mean_q_clean | Float | Mean quality score of clean forward reads |
| r1_mean_q_raw | Float | Mean quality score of raw forward reads |
| r2_mean_q_clean | Float | Mean quality score of clean reverse reads |
| r2_mean_q_raw | Float | Mean quality score of raw reverse reads |
| fastq_scan_version | String | Version of fastq-scan software used |
| gambit_closest_genomes | File | CSV file listing genomes in the GAMBIT database that are most similar to the query assembly |
| gambit_db_version | String | Version of GAMBIT used |
| gambit_docker | String | GAMBIT docker file used |
| gambit_predicted_taxon | String | Taxon predicted by GAMBIT |
| gambit_predicted_taxon_rank | String | Taxon rank of GAMBIT taxon prediction |
| gambit_report | File | GAMBIT report in a machine-readable format |
| gambit_version | String | Version of GAMBIT software used |
| assembly_length | Int | Length of assembly (total contig length) as determined by QUAST |
| n50_value | Int | N50 of assembly calculated by QUAST |
| number_contigs | Int | Total number of contigs in assembly |
| quast_report | File | TSV report from QUAST |
| quast_version | String | Software version of QUAST used |
| rasusa_version | String | Version of rasusa used |
| read1_subsampled | File | Subsampled read1 file |
| read2_subsampled | File | Subsampled read2 file |
| bbduk_docker | String | BBDuk docker image used  |
| fastp_version | String | Version of fastp software used |
| read1_clean | File | Clean forward reads file |
| read2_clean | File | Clean reverse reads file |
| num_reads_clean_pairs | String | Number of read pairs after cleaning |
| num_reads_clean1 | Int | Number of forward reads after cleaning |
| num_reads_clean2 | Int | Number of reverse reads after cleaning |
| num_reads_raw_pairs | String | Number of input read pairs |
| num_reads_raw1 | Int | Number of input forward reads |
| num_reads_raw2 | Int | Number of input reverse reads |
| trimmomatic_version | String | Version of trimmomatic used |
| clean_read_screen | String | PASS or FAIL result from clean read screening; FAIL accompanied by the reason for failure |
| raw_read_screen | String | PASS or FAIL result from raw read screening; FAIL accompanied by thereason for failure |
| assembly_fasta | File | <https://github.com/tseemann/shovill#contigsfa> |
| contigs_fastg | File | Assembly graph if megahit used for genome assembly |
| contigs_gfa | File | Assembly graph if spades used for genome assembly |
| contigs_lastgraph | File | Assembly graph if velvet used for genome assembly |
| shovill_pe_version | String | Shovill version used |
| theiaeuk_snippy_variants_bam | File | BAM file produced by the snippy module |
| theiaeuk_snippy_variants_gene_query_results | File | File containing all lines from variants file matching gene query terms |
| theiaeuk_snippy_variants_hits | String | String of all variant file entries matching gene query term |
| theiaeuk_snippy_variants_outdir_tarball | File | Tar compressed file containing full snippy output directory |
| theiaeuk_snippy_variants_query | String | The gene query term(s) used to search variant |
| theiaeuk_snippy_variants_query_check | String | Were the gene query terms present in the refence annotated genome file |
| theiaeuk_snippy_variants_reference_genome | File | The reference genome used in the alignment and variant calling |
| theiaeuk_snippy_variants_results | File | The variants file produced by snippy |
| theiaeuk_snippy_variants_summary | File | A file summarizing the variants detected by snippy |
| theiaeuk_snippy_variants_version | String | The version of the snippy_variants module being used |
| seq_platform | String | Sequencing platform inout by the user |
| theiaeuk_illumina_pe_analysis_date | String | Date of TheiaProk workflow execution |
| theiaeuk_illumina_pe_version | String | TheiaProk workflow version used |
