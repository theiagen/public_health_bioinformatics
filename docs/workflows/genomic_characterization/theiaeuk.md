# TheiaEuk Workflow Series

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibliity** | **Workflow Level** |
|---|---|---|---|---|
| [Genomic Characterization](../../workflows_overview/workflows_type.md/#genomic-characterization) | [Mycotics](../../workflows_overview/workflows_kingdom.md/#mycotics) | PHB v3.0.1 | Yes | Sample-level |

## TheiaEuk Workflows

**The TheiaEuk_Illumina_PE workflow is for the assembly, quality assessment, and characterization of fungal genomes.** It is designed to accept Illumina paired-end sequencing data as the primary input. **It is currently intended only for ==haploid== fungal genomes like _Candidozyma auris_.** Analyzing diploid genomes using TheiaEuk should be attempted only with expert attention to the resulting genome quality.

All input reads are processed through "core tasks" in each workflow. The core tasks include raw read quality assessment, read cleaning (quality trimming and adapter removal), de novo assembly, assembly quality assessment, and species taxon identification. For some taxa identified, taxa-specific sub-workflows will be automatically activated, undertaking additional taxa-specific characterization steps, including clade-typing and/or antifungal resistance detection.

!!! caption "TheiaEuk Workflow Diagram"
    ![TheiaEuk Workflow Diagram](../../assets/figures/TheiaEuk_Illumina_PHB_202557.png){width=75%}

### Inputs

!!! info "Input read data"

    The TheiaEuk_Illumina_PE workflow takes in Illumina paired-end read data. Read file names should end with `.fastq` or `.fq`, with the optional addition of `.gz`. When possible, Theiagen recommends zipping files with [gzip](https://www.gnu.org/software/gzip/) prior to Terra upload to minimize data upload time.

    By default, the workflow anticipates 2 x 150bp reads (i.e. the input reads were generated using a 300-cycle sequencing kit). Modifications to the optional parameter for `trim_minlen` may be required to accommodate shorter read data, such as the 2 x 75bp reads generated using a 150-cycle sequencing kit.

<div class="searchable-table" markdown="1">

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
| digger_denovo | **assembler** | String | Assembler to use (spades, skesa, megahit) | spades | Optional |
| digger_denovo | **assembly_options** | String | String | Assembler-specific options that you might choose for the selected assembler | | Optional |
| digger_denovo | **bwa_cpu** | Int | Number of CPU cores for BWA alignment | 6 | Optional |
| digger_denovo | **bwa_disk_size** | Int | Disk space in GB for BWA alignment | 100 | Optional |
| digger_denovo | **bwa_docker** | String | Docker image for BWA alignment | us-docker.pkg.dev/general-theiagen/staphb/ivar:1.3.1-titan | Optional |
| digger_denovo | **bwa_memory** | Int | Memory in GB for BWA alignment | 16 | Optional |
| digger_denovo | **filter_contigs_cpu** | Int | Number of CPU cores for contig filtering | 1 | Optional |
| digger_denovo | **filter_contigs_disk_size** | Int | Disk space in GB for contig filtering | 50 | Optional |
| digger_denovo | **filter_contigs_docker** | String | Docker image for contig filtering | us-docker.pkg.dev/general-theiagen/theiagen/shovilter:0.2 | Optional |
| digger_denovo | **filter_contigs_memory** | Int | Memory in GB for contig filtering | 8 | Optional |
| digger_denovo | **filter_contigs_min_coverage** | Float | Minimum coverage threshold for contig filtering | 5.0 | Optional |
| digger_denovo | **filter_contigs_min_length** | Int | Minimum length threshold for contig filtering | 200 | Optional |
| digger_denovo | **filter_contigs_skip_coverage_filter** | Boolean | Skip filtering contigs based on coverage | false | Optional |
| digger_denovo | **filter_contigs_skip_homopolymer_filter** | Boolean | Skip filtering contigs containing homopolymers | false | Optional |
| digger_denovo | **filter_contigs_skip_length_filter** | Boolean | Skip filtering contigs based on length | false | Optional |
| digger_denovo | **kmers** | String | K-mer sizes for assembly (comma-separated) | | Optional |
| digger_denovo | **megahit_cpu** | Int | Number of CPU cores for MEGAHIT assembler | 4 | Optional |
| digger_denovo | **megahit_disk_size** | Int | Disk space in GB for MEGAHIT assembler | 100 | Optional |
| digger_denovo | **megahit_docker** | String | Docker image for MEGAHIT assembler | us-docker.pkg.dev/general-theiagen/theiagen/megahit:1.2.9 | Optional |
| digger_denovo | **megahit_memory** | Int | Memory in GB for MEGAHIT assembler | 16 | Optional |
| digger_denovo | **min_contig_length** | Int | Minimum contig length to retain in final assembly | 200 | Optional |
| digger_denovo | **pilon_cpu** | Int | Number of CPU cores for Pilon polishing | 8 | Optional |
| digger_denovo | **pilon_disk_size** | Int | Disk space in GB for Pilon polishing | 100 | Optional |
| digger_denovo | **pilon_docker** | String | Docker image for Pilon polishing | us-docker.pkg.dev/general-theiagen/biocontainers/pilon:1.24--hdfd78af_0 | Optional |
| digger_denovo | **pilon_fix** | String | Potential issues with assembly to try and automatically fix (snps, indels, gaps, local, all, bases, none) | "bases" | Optional |
| digger_denovo | **pilon_memory** | Int | Memory in GB for Pilon polishing | 32 | Optional |
| digger_denovo | **pilon_min_base_quality** | Int | Minimum base quality to keep | 3 | Optional |
| digger_denovo | **pilon_min_depth** | Float | Minimum coverage threshold for variant calling: when set to a value ≥1, it requires that absolute depth of coverage; when set to a fraction <1, it requires coverage at least that fraction of the mean coverage for the region | 0.25 | Optional |
| digger_denovo | **pilon_min_mapping_quality** | Int | Minimum mapping quality for a read to count in pileups | 60 |
| digger_denovo | **run_filter_contigs** | Boolean | Whether to run contig filtering step | true | Optional |
| digger_denovo | **skesa_cpu** | Int | Number of CPU cores for SKESA assembler | 4 | Optional |
| digger_denovo | **skesa_disk_size** | Int | Disk space in GB for SKESA assembler | 50 | Optional |
| digger_denovo | **skesa_docker** | String | Docker image for SKESA assembler | us-docker.pkg.dev/general-theiagen/staphb/skesa:2.4.0 | Optional |
| digger_denovo | **skesa_memory** | Int | Memory in GB for SKESA assembler | 4 | Optional |
| digger_denovo | **spades_cpu** | Int | Number of CPU cores for SPAdes assembler | 16 | Optional |
| digger_denovo | **spades_disk_size** | Int | Disk space in GB for SPAdes assembler | 100 | Optional |
| digger_denovo | **spades_docker** | String | Docker image for SPAdes assembler | us-docker.pkg.dev/general-theiagen/staphb/spades:4.1.0 | Optional |
| digger_denovo | **spades_memory** | Int | Memory in GB for SPAdes assembler | 32 | Optional |
| digger_denovo | **spades_type** | String | SPAdes assembly mode (isolate, meta, rna, etc.), more can be found [here](https://ablab.github.io/spades/running.html) | isolate | Optional |
| digger_denovo | **use_pilon** | Boolean | Whether to run Pilon polishing after assembly | false | Optional |
| gambit | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| gambit | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/gambit:1.0.0 | Optional |
| merlin_magic | **agrvate_docker_image** | String | Internal component, do not modify | "us-docker.pkg.dev/general-theiagen/biocontainers/agrvate:1.0.2--hdfd78af_0" | Do Not Modify, Optional |
| merlin_magic | **assembly_only** | Boolean | Internal component, do not modify |  | Do Not Modify, Optional |
| merlin_magic | **amr_search_cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| merlin_magic | **amr_search_disk_size** | Int | Amount of storage (in GB) to allocate to the task | 50 | Optional |
| merlin_magic | **amr_search_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/amrsearch:0.2.1 | Optional |
| merlin_magic | **amr_search_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| merlin_magic | **call_poppunk** | Boolean | Internal component, do not modify | TRUE | Do Not Modify, Optional |
| merlin_magic | **call_shigeifinder_reads_input** | Boolean | Internal component, do not modify | FALSE | Do Not Modify, Optional |
| merlin_magic | **emmtypingtool_docker_image** | String | Internal component, do not modify | us-docker.pkg.dev/general-theiagen/staphb/emmtypingtool:0.0.1 | Do Not Modify, Optional |
| merlin_magic | **hicap_docker_image** | String | Internal component, do not modify | us-docker.pkg.dev/general-theiagen/biocontainers/hicap:1.0.3--py_0 | Do Not Modify, Optional |
| merlin_magic | **ont_data** | Boolean | Internal component, do not modify |  | Do Not Modify, Optional |
| merlin_magic | **paired_end** | Boolean | Internal component, do not modify |  | Do Not Modify, Optional |
| merlin_magic | **pasty_docker_image** | String | Internal component, do not modify | us-docker.pkg.dev/general-theiagen/staphb/pasty:1.0.3 | Do Not Modify, Optional |
| merlin_magic | **pasty_min_coverage** | Int | Internal component, do not modify | 95 | Do Not Modify, Optional |
| merlin_magic | **pasty_min_percent_identity** | Int | Internal component, do not modify | 95 | Do Not Modify, Optional |
| merlin_magic | **run_amr_search** | Boolean | If set to true AMR_Search workflow will be run if species is part of supported taxon, see AMR_Search docs. | False | Optional |
| merlin_magic | **shigatyper_docker_image** | String | Internal component, do not modify | us-docker.pkg.dev/general-theiagen/staphb/shigatyper:2.0.5 | Do Not Modify, Optional |
| merlin_magic | **shigeifinder_docker_image** | String | Internal component, do not modify | us-docker.pkg.dev/general-theiagen/staphb/shigeifinder:1.3.5 | Do Not Modify, Optional |
| merlin_magic | **snippy_query_gene** | String | Internal component, do not modify |  | Do Not Modify, Optional |
| merlin_magic | **srst2_gene_max_mismatch** | Int | Internal component, do not modify | 2000 | Do Not Modify, Optional |
| merlin_magic | **srst2_max_divergence** | Int | Internal component, do not modify | 20 | Do Not Modify, Optional |
| merlin_magic | **srst2_min_cov** | Int | Internal component, do not modify | 80 | Do Not Modify, Optional |
| merlin_magic | **srst2_min_depth** | Int | Internal component, do not modify | 5 | Do Not Modify, Optional |
| merlin_magic | **srst2_min_edge_depth** | Int | Internal component, do not modify | 2 | Do Not Modify, Optional |
| merlin_magic | **staphopia_sccmec_docker_image** | String | Internal component, do not modify | us-docker.pkg.dev/general-theiagen/biocontainers/staphopia-sccmec:1.0.0--hdfd78af_0 | Do Not Modify, Optional |
| merlin_magic | **tbp_parser_config** | File | Internal component, do not modify |  | Do Not Modify, Optional |
| merlin_magic | **tbp_parser_min_percent_coverage** | Float | Internal component, do not modify | 100.0 | Do Not Modify, Optional |
| merlin_magic | **tbp_parser_debug** | Boolean | Internal component, do not modify | FALSE | Do Not Modify, Optional |
| merlin_magic | **tbp_parser_docker_image** | String | Internal component, do not modify | us-docker.pkg.dev/general-theiagen/theiagen/tbp-parser:2.4.5 | Do Not Modify, Optional |
| merlin_magic | **tbp_parser_min_depth** | Int | Internal component, do not modify | 10 | Do Not Modify, Optional |
| merlin_magic | **tbp_parser_operator** | String | Internal component, do not modify | "Operator not provided" | Do Not Modify, Optional |
| merlin_magic | **tbp_parser_output_seq_method_type** | String | Internal component, do not modify | "WGS" | Do Not Modify, Optional |
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
| merlin_magic | **virulencefinder_database** | String | Internal component, do not modify | "virulence_ecoli" | Do Not Modify, Optional |
| merlin_magic | **virulencefinder_docker_image** | String | Internal component, do not modify | us-docker.pkg.dev/general-theiagen/staphb/virulencefinder:2.0.4 | Do Not Modify, Optional |
| merlin_magic | **virulencefinder_min_percent_coverage** | Float | Internal component, do not modify |  | Do Not Modify, Optional |
| merlin_magic | **virulencefinder_min_percent_identity** | Float | Internal component, do not modify |  | Do Not Modify, Optional |
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
| rasusa_task | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/rasusa:2.1.0 | Optional |
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
| theiaeuk_pe | **busco_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| theiaeuk_pe | **call_rasusa** | Boolean | If true, launch rasusa task to subsample raw reads to read depth of 150X | TRUE | Optional |
| theiaeuk_pe | **gambit_db_genomes** | File | User-provided database of assembled query genomes; requires complementary signatures file. If not provided, uses default database, "/gambit-db" | gs://gambit-databases-rp/fungal-version/1.0.0/gambit-fungal-metadata-1.0.0-20241213.gdb | Optional |
| theiaeuk_pe | **gambit_db_signatures** | File | User-provided signatures file; requires complementary genomes file. If not specified, the file from the docker container will be used.  | gs://gambit-databases-rp/fungal-version/1.0.0/gambit-fungal-signatures-1.0.0-20241213.gs | Optional |
| theiaeuk_pe | **genome_length** | Int | User-specified expected genome size to be used in genome statistics calculations |  | Optional |
| theiaeuk_pe | **max_genome_size** | Int | Maximum genome size able to pass read screening | 50000000 | Optional |
| theiaeuk_pe | **min_basepairs** | Int | Minimum number of base pairs able to pass read screening | 2241820 | Optional |
| theiaeuk_pe | **min_coverage** | Int | Minimum genome coverage able to pass read screening | 10 | Optional |
| theiaeuk_pe | **min_genome_size** | Int | Minimum genome size able to pass read screening | 100000 | Optional |
| theiaeuk_pe | **min_proportion** | Int | Minimum proportion of total reads in each read file to pass read screening | 50 | Optional |
| theiaeuk_pe | **min_reads** | Int | Minimum number of reads to pass read screening | 10000 | Optional |
| theiaeuk_pe | **skip_screen** | Boolean | Option to skip the read screening prior to analysis; if setting to true, please provide a value for the theiaeuk_pe `genome_length` optional input, OR set `call_rasusa` to false. Otherwise RASUSA will attempt to downsample to an expected genome size of 0 bp, and the workflow will fail. | FALSE | Optional |
| theiaeuk_pe | **subsample_coverage** | Float | Read depth for RASUSA task to subsample reads to | 150 | Optional |
| version_capture | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional |

</div>

### Workflow Tasks

All input reads are processed through "core tasks" in the TheiaEuk workflows. These undertake read trimming and assembly appropriate to the input data type, currently only Illumina paired-end data. TheiaEuk workflow subsequently launch default genome characterization modules for quality assessment, and additional taxa-specific characterization steps. When setting up the workflow, users may choose to use "optional tasks" or alternatives to tasks run in the workflow by default.

#### Core tasks

!!! tip ""
    These tasks are performed regardless of organism. They perform read trimming and various quality control steps.

??? task "`versioning`: Version capture for TheiaEuk"

    The `versioning` task captures the workflow version from the GitHub (code repository) version.
        
    !!! techdetails "Version Capture Technical details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_versioning.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/task_versioning.wdl) |

??? task "`screen`: Total Raw Read Quantification and Genome Size Estimation (optional, on by default)"

    The [`screen`](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/comparisons/task_screen.wdl) task ensures the quantity of sequence data is sufficient to undertake genomic analysis. It uses [`fastq-scan`](https://github.com/rpetit3/fastq-scan) and bash commands for quantification of reads and base pairs, and [mash](https://mash.readthedocs.io/en/latest/index.html) sketching to estimate the genome size and its coverage. At each step, the results are assessed relative to pass/fail criteria and thresholds that may be defined by optional user inputs. Samples are run through all threshold checks, regardless of failures, and the workflow will terminate after the `screen` task if any thresholds are not met:


    1. Total number of reads: A sample will fail the read screening task if its total number of reads is less than or equal to `min_reads`.
    2. The proportion of basepairs reads in the forward and reverse read files: A sample will fail the read screening if fewer than `min_proportion` basepairs are in either the reads1 or read2 files.
    3. Number of basepairs: A sample will fail the read screening if there are fewer than `min_basepairs` basepairs
    4. Estimated genome size:  A sample will fail the read screening if the estimated genome size is smaller than `min_genome_size` or bigger than `max_genome_size`.
    5. Estimated genome coverage: A sample will fail the read screening if the estimated genome coverage is less than the `min_coverage`.

    Read screening is undertaken on both the raw and cleaned reads. The task may be skipped by setting the `skip_screen` variable to true.

    Default values vary between the PE and SE workflow. The rationale for these default values can be found below.
    
    | Variable  | Rationale |
    | --- | --- |
    | `skip_screen` | Prevent the read screen from running. If you set this value to true, please provide a value for the theiaeuk_pe `genome_length` optional input, OR set the theiaeuk_pe `call_rasusa` optional input to false. Otherwise RASUSA will attempt to downsample to an expected genome size of 0 bp, and the workflow will fail. |
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
        | Task | [task_screen.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/comparisons/task_screen.wdl)  |

??? task "`Rasusa`: Read subsampling (optional, on by default)"

    The Rasusa task performs subsampling of the raw reads. By default, this task will subsample reads to a depth of 150X using the estimated genome length produced during the preceding raw read screen. The user can prevent the task from being launched by setting the `call_rasusa`variable to false. 

    The user can also provide an estimated genome length for the task to use for subsampling using the `genome_size` variable. In addition, the read depth can be modified using the `subsample_coverage` variable.
        
    !!! techdetails "Rasusa Technical Details"

        |  | Links |
        | --- | --- |
        | Task | [task_rasusa.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/task_rasusa.wdl) |
        | Software Source Code | [Rasusa on GitHub](https://github.com/mbhall88/rasusa) |
        | Software Documentation | [Rasusa on GitHub](https://github.com/mbhall88/rasusa) |
        | Original Publication(s) | [Rasusa: Randomly subsample sequencing reads to a specified coverage](https://doi.org/10.21105/joss.03941) |

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
        | Sub-workflow | [wf_read_QC_trim_pe.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_read_QC_trim_pe.wdl) |
        | Tasks | [task_fastp.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_fastp.wdl)<br>[task_trimmomatic.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_trimmomatic.wdl)<br>[task_bbduk.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_bbduk.wdl)<br>[task_fastq_scan.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_fastq_scan.wdl)<br>[task_midas.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/contamination/task_midas.wdl)<br>[task_kraken2.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/contamination/task_kraken2.wdl)|
        | Software Source Code | [fastp](https://github.com/OpenGene/fastp); [Trimmomatic](https://github.com/usadellab/Trimmomatic); [fastq-scan](https://github.com/rpetit3/fastq-scan); [MIDAS](https://github.com/snayfach/MIDAS); [Kraken2](https://github.com/DerrickWood/kraken2)|
        | Software Documentation | [fastp](https://github.com/OpenGene/fastp); [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic); [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/); [fastq-scan](https://github.com/rpetit3/fastq-scan); [MIDAS](https://github.com/snayfach/MIDAS); [Kraken2](https://github.com/DerrickWood/kraken2/wiki) |
        | Original Publication(s) | [Trimmomatic: a flexible trimmer for Illumina sequence data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4103590/)<br>[fastp: an ultra-fast all-in-one FASTQ preprocessor](https://academic.oup.com/bioinformatics/article/34/17/i884/5093234?login=false)<br>[An integrated metagenomics pipeline for strain profiling reveals novel patterns of bacterial transmission and biogeography](https://pubmed.ncbi.nlm.nih.gov/27803195/)<br>[Improved metagenomic analysis with Kraken 2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0) |

#### Assembly tasks

!!! tip ""
    These tasks assemble the reads into a _de novo_ assembly and assess the quality of the assembly.

??? task "`digger_denovo`: _De novo_ Assembly"

    De Novo assembly will be undertaken only for samples that have sufficient read quantity and quality, as determined by the `screen` task assessment of clean reads. 

    In TheiaEuk, assembly is performed using the [digger_denovo](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_digger_denovo.wdl), which is a hat tip to [Shovill](https://github.com/tseemann/shovill) pipeline. This undertakes the assembly with one of three assemblers [SKESA](https://github.com/ncbi/SKESA) (default), [SPAdes](https://github.com/ablab/spades, [Megahit](https://github.com/voutcn/megahit)), but also performs a number of post processing steps for assembly polishing and contig filtering. Pilon can optionally be run if `use_pilon` is set to true. On defualt, the contig filtering task is set to run, which will remove any homopolymers, contigs below a specificied length, and contigs with coverage below a specified minimum coverage. This can be turned off by setting `run_filter_contigs` to `false`. 
    
    ??? toggle "What is _de novo_  assembly?"
        _De novo_  assembly is the process or product of attempting to reconstruct a genome from scratch (without prior knowledge of the genome) using sequence reads. Assembly of fungal genomes from short-reads will produce multiple contigs per chromosome rather than a single contiguous sequence for each chromosome.
        
    !!! techdetails "Digger-Denovo Technical Details"
        |  | Links |
        | --- | --- |
        | TheiaEuk WDL SubWorkflow | [wf_digger_denovo.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_digger_denovo.wd) |
        | Software Source Code | [digger_denovo on GitHub](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_digger_denovo.wd) |

??? task "`QUAST`: Assembly Quality Assessment"

    [`QUAST`](https://github.com/ablab/quast) (**QU**ality **AS**sessment **T**ool) evaluates genome assemblies by computing several metrics that describe the assembly quality, including the total number of bases in the assembly, the length of the largest contig in the assembly, and the assembly percentage GC content.

    !!! techdetails "QUAST Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_quast.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_quast.wdl) |
        | Software Source Code | [QUAST on GitHub](https://github.com/ablab/quast) |
        | Software Documentation | https://quast.sourceforge.net/docs/manual.html |
        | Orginal publication | [QUAST: quality assessment tool for genome assemblies](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3624806/) |

??? task "`CG-Pipeline`: Assessment of Read Quality, and Estimation of Genome Coverage"

    The`cg_pipeline` task generates metrics about read quality and estimates the coverage of the genome using the "run_assembly_readMetrics.pl" script from [CG-Pipeline](https://github.com/lskatz/CG-Pipeline/). The genome coverage estimates are calculated using both using raw and cleaned reads, using either a user-provided `genome_size` or the estimated genome length generated by QUAST.

    !!! techdetails "CG-Pipeline Technical Details"
        The `cg_pipeline` task is run twice in TheiaEuk, once with raw reads, and once with clean reads.
        
        |  | Links |
        | --- | --- |
        | Task | [task_cg_pipeline.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_cg_pipeline.wdl) |
        | Software Source Code | [CG-Pipeline on GitHub](https://github.com/lskatz/CG-Pipeline/) |
        | Software Documentation | [CG-Pipeline on GitHub](https://github.com/lskatz/CG-Pipeline/) |
        | Original Publication(s) | [A computational genomics pipeline for prokaryotic sequencing projects](https://academic.oup.com/bioinformatics/article/26/15/1819/188418) |

#### Organism-agnostic characterization

!!! tip ""
    These tasks are performed regardless of the organism and provide quality control and taxonomic assignment.

??? task "`GAMBIT`: **Taxon Assignment**"

    [`GAMBIT`](https://github.com/jlumpe/gambit) determines the taxon of the genome assembly using a k-mer based approach to match the assembly sequence to the closest complete genome in a database, thereby predicting its identity. Sometimes, GAMBIT can confidently designate the organism to the species level. Other times, it is more conservative and assigns it to a higher taxonomic rank.

    For additional details regarding the GAMBIT tool and a list of available GAMBIT databases for analysis, please consult the [GAMBIT](../../guides/gambit.md) tool documentation.

    !!! techdetails "GAMBIT Technical Details"

        |  | Links |
        | --- | --- |
        | Task | [task_gambit.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/task_gambit.wdl) |
        | Software Source Code | [GAMBIT on GitHub](https://github.com/jlumpe/gambit) |
        | Software Documentation | [GAMBIT ReadTheDocs](https://gambit-genomics.readthedocs.io/en/latest/) |
        | Original Publication(s) | [GAMBIT (Genomic Approximation Method for Bacterial Identification and Tracking): A methodology to rapidly leverage whole genome sequencing of bacterial isolates for clinical identification](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0277575) |

??? task "`BUSCO`: Assembly Quality Assessment"

    BUSCO (**B**enchmarking **U**niversal **S**ingle-**C**opy **O**rthologue) attempts to quantify the completeness and contamination of an assembly to generate quality assessment metrics. It uses taxa-specific databases containing genes that are all expected to occur in the given taxa, each in a single copy. BUSCO examines the presence or absence of these genes, whether they are fragmented, and whether they are duplicated (suggestive that additional copies came from contaminants).

    **BUSCO notation** 
    
    Here is an example of BUSCO notation: `C:99.1%[S:98.9%,D:0.2%],F:0.0%,M:0.9%,n:440`. There are several abbreviations used in this output:
    
    - Complete (C) - genes are considered "complete" when their lengths are within two standard deviations of the BUSCO group mean length.
    - Single-copy (S) - genes that are complete and have only one copy.
    - Duplicated (D) - genes that are complete and have more than one copy.
    - Fragmented (F) - genes that are only partially recovered.
    - Missing (M) - genes that were not recovered at all.
    - Number of genes examined (n) - the number of genes examined.
    
    A high equity assembly will use the appropriate database for the taxa, have high complete (C) and single-copy (S) percentages, and low duplicated (D), fragmented (F) and missing (M) percentages. 
  
    !!! techdetails "BUSCO Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_busco.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/advanced_metrics/task_busco.wdl) |
        | Software Source Code | [BUSCO on GitLab](https://gitlab.com/ezlab/busco) |
        | Software Documentation | https://busco.ezlab.org/ |
        | Orginal publication | [BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs](https://academic.oup.com/bioinformatics/article/31/19/3210/211866) |

??? task "`qc_check`: Check QC Metrics Against User-Defined Thresholds (optional)"

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
        | Task | [task_qc_check_phb.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/comparisons/task_qc_check_phb.wdl) |

#### Organism-specific characterization

!!! tip ""
    The TheiaEuk workflow automatically activates taxa-specific tasks after identification of the relevant taxa using `GAMBIT`. Many of these taxa-specific tasks do not require any additional inputs from the user.

??? toggle "_Candidozyma auris_ (also known as _Candida auris_)"
    Two tools are deployed when _Candidozyma auris_/_Candida auris_ is  identified.

    ??? task "Cladetyping: clade determination"
        GAMBIT is used to determine the clade of the specimen by comparing the sequence to five clade-specific reference files. The output of the clade typing task will be used to specify the reference genome for the antifungal resistance detection tool.

        ??? toggle "Default reference genomes used for clade typing and antimicrobial resistance gene detection of _C. auris_"
            | Clade | Genome Accession | Assembly Name | Strain | NCBI Submitter | Included mutations in AMR genes (not comprehensive) |
            | --- | --- | --- | --- | --- | --- |
            | _Candidozyma auris_ Clade I | GCA_002759435.2 | Cand_auris_B8441_V2 | B8441 | Centers for Disease Control and Prevention |  |
            | _Candidozyma auris_ Clade II | GCA_003013715.2 | ASM301371v2 | B11220 | Centers for Disease Control and Prevention |  |
            | _Candidozyma auris_ Clade III | GCA_002775015.1 | Cand_auris_B11221_V1 | B11221 | Centers for Disease Control and Prevention | _ERG11_ V125A/F126L |
            | _Candidozyma auris_ Clade IV | GCA_003014415.1 | Cand_auris_B11243 | B11243 | Centers for Disease Control and Prevention | _ERG11_ Y132F |
            | _Candidozyma auris_ Clade V | GCA_016809505.1 | ASM1680950v1 | IFRC2087 | Centers for Disease Control and Prevention |  |

        !!! techdetails "Cladetyping Technical Details"
            |  | Links |
            | --- | --- |
            | Task | [task_cauris_cladetyping.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/candida/task_cauris_cladetyper.wdl) |
            | Software Source Code | [GAMBIT on GitHub](https://github.com/jlumpe/gambit) |
            | Software Documentation | [GAMBIT Overview](https://theiagen.notion.site/GAMBIT-7c1376b861d0486abfbc316480046bdc?pvs=4)
            | Original Publication(s) | [GAMBIT (Genomic Approximation Method for Bacterial Identification and Tracking): A methodology to rapidly leverage whole genome sequencing of bacterial isolates for clinical identification](https://doi.org/10.1371/journal.pone.0277575)<br> [TheiaEuk: a species-agnostic bioinformatics workflow for fungal genomic characterization](https://doi.org/10.3389/fpubh.2023.1198213) |

    ??? task "Snippy Variants: antifungal resistance detection"
        To detect mutations that may confer antifungal resistance, `Snippy` is used to find all variants relative to the clade-specific reference, then these variants are queried for product names associated with resistance.
    
        The genes in which there are known resistance-conferring mutations for this pathogen are:

        - FKS1
        - ERG11 (lanosterol 14-alpha demethylase)
        - FUR1 (uracil phosphoribosyltransferase)

        We query `Snippy` results to see if any mutations were identified in those genes. By default, we automatically check for the following loci (which can be overwritten by the user). You will find the mutations next to the locus tag in the `theiaeuk_snippy_variants_hits` column corresponding gene name (see below):

        | **TheiaEuk Search Term** | **Corresponding Gene Name** |
        |---|---|
        | B9J08_005340 | ERG6 |
        | B9J08_000401 | FLO8 |
        | B9J08_005343 | Hypothetical protein (PSK74852) |
        | B9J08_003102 | MEC3 |
        | B9J08_003737 | ERG3 |
        | lanosterol.14-alpha.demethylase | ERG11 |
        | uracil.phosphoribosyltransferase | FUR1 |
        | FKS1 | FKS1 |    

        For example, one sample may have the following output for the `theiaeuk_snippy_variants_hits` column:

        ```plaintext
        lanosterol.14-alpha.demethylase: lanosterol 14-alpha demethylase (missense_variant c.428A>G p.Lys143Arg; C:266 T:0),B9J08_000401: hypothetical protein (stop_gained c.424C>T p.Gln142*; A:70 G:0)
        ```

        Based on this, we can tell that ERG11 has a missense variant at position 143 (Lysine to Arginine) and B9J08_000401 (which is FLO8) has a stop-gained variant at position 142 (Glutamine to Stop).

        ??? toggle "Known resistance-conferring mutations for _Candidozyma auris_"
            Mutations in these genes that are known to confer resistance are shown below

            | **Organism** | **Found in** | **Gene name** | **Gene locus** | **AA mutation** | **Drug** | **Reference** |
            | --- | --- | --- | --- | --- | --- | --- |
            | **Candidozyma auris** | **Human** | **ERG11** |  | **Y132F** | **Fluconazole** | [Simultaneous Emergence of Multidrug-Resistant _Candida auris_ on 3 Continents Confirmed by Whole-Genome Sequencing and Epidemiological Analyses](https://academic.oup.com/cid/article/64/2/134/2706620/Simultaneous-Emergence-of-Multidrug-Resistant) |
            | **Candidozyma auris** | **Human** | **ERG11** |  | **K143R** | **Fluconazole** | [Simultaneous Emergence of Multidrug-Resistant _Candida auris_ on 3 Continents Confirmed by Whole-Genome Sequencing and Epidemiological Analyses](https://academic.oup.com/cid/article/64/2/134/2706620/Simultaneous-Emergence-of-Multidrug-Resistant) |
            | **Candidozyma auris** | **Human** | **ERG11** |  | **F126T** | **Fluconazole** | [Simultaneous Emergence of Multidrug-Resistant _Candida auris_ on 3 Continents Confirmed by Whole-Genome Sequencing and Epidemiological Analyses](https://academic.oup.com/cid/article/64/2/134/2706620/Simultaneous-Emergence-of-Multidrug-Resistant) |
            | **Candidozyma auris** | **Human** | **FKS1** |  | **S639P** | **Micafungin**  | [Activity of CD101, a long-acting echinocandin, against clinical isolates of Candida auris](https://www.sciencedirect.com/science/article/pii/S0732889317303498) |
            | **Candidozyma auris** | **Human** | **FKS1** |  | **S639P** | **Caspofungin** | [Activity of CD101, a long-acting echinocandin, against clinical isolates of Candida auris](https://www.sciencedirect.com/science/article/pii/S0732889317303498) |
            | **Candidozyma auris** | **Human** | **FKS1** |  | **S639P** | **Anidulafungin** | [Activity of CD101, a long-acting echinocandin, against clinical isolates of Candida auris](https://www.sciencedirect.com/science/article/pii/S0732889317303498) |
            | **Candidozyma auris** | **Human** | **FKS1** |  | **S639F** | **Micafungin** | [A multicentre study of antifungal susceptibility patterns among 350 _Candida auris_ isolates (2009–17) in India: role of the ERG11 and FKS1 genes in azole and echinocandin resistance](https://academic.oup.com/jac/advance-article/doi/10.1093/jac/dkx480/4794718) |
            | **Candidozyma auris** | **Human** | **FKS1** |  | **S639F** | **Caspofungin** | [A multicentre study of antifungal susceptibility patterns among 350 _Candida auris_ isolates (2009–17) in India: role of the ERG11 and FKS1 genes in azole and echinocandin resistance](https://academic.oup.com/jac/advance-article/doi/10.1093/jac/dkx480/4794718) |
            | **Candidozyma auris** | **Human** | **FKS1** |  | **S639F** | **Anidulafungin** | [A multicentre study of antifungal susceptibility patterns among 350 _Candida auris_ isolates (2009–17) in India: role of the ERG11 and FKS1 genes in azole and echinocandin resistance](https://academic.oup.com/jac/advance-article/doi/10.1093/jac/dkx480/4794718) |
            | **Candidozyma auris** | **Human** | **FUR1** | **CAMJ_004922** | **F211I** | **5-flucytosine** | [Genomic epidemiology of the UK outbreak of the emerging human fungal pathogen Candida auris](https://doi.org/10.1038/s41426-018-0045-x) |

        !!! techdetails "Snippy Variants Technical Details"
            |  | Links |
            | --- | --- |
            | Task | [task_snippy_variants.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snippy_variants.wdl)<br>[task_snippy_gene_query.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snippy_gene_query.wdl) |
            | Software Source Code | [Snippy on GitHub](https://github.com/tseemann/snippy) |
            | Software Documentation | [Snippy on GitHub](https://github.com/tseemann/snippy) |
            
??? toggle "_Candida albicans_"
    When this species is detected by the taxon ID tool, an antifungal resistance detection task is deployed.

    ??? task "Snippy Variants: antifungal resistance detection"
        To detect mutations that may confer antifungal resistance, `Snippy` is used to find all variants relative to the clade-specific reference, and these variants are queried for product names associated with resistance.

        The genes in which there are known resistance-conferring mutations for this pathogen are:

        - ERG11
        - GCS1 (FKS1)
        - FUR1
        - RTA2

        We query `Snippy` results to see if any mutations were identified in those genes. By default, we automatically check for the following loci (which can be overwritten by the user). You will find the mutations next to the locus tag in the `theiaeuk_snippy_variants_hits` column corresponding gene name (see below):

        | **TheiaEuk Search Term** | **Corresponding Gene Name** |
        |---|---|
        | ERG11 | ERG11 |
        | GCS1 | FKS1 |
        | FUR1 | FUR1 |
        | RTA2 | RTA2 |

        !!! techdetails "Snippy Variants Technical Details"
            |  | Links |
            | --- | --- |
            | Task | [task_snippy_variants.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snippy_variants.wdl)<br>[task_snippy_gene_query.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snippy_gene_query.wdl) |
            | Software Source Code | [Snippy on GitHub](https://github.com/tseemann/snippy) |
            | Software Documentation | [Snippy on GitHub](https://github.com/tseemann/snippy) |

??? toggle "_Aspergillus fumigatus_"
    When this species is detected by the taxon ID tool an antifungal resistance detection task is deployed.

    ??? task "Snippy Variants: antifungal resistance detection"
        To detect mutations that may confer antifungal resistance, `Snippy` is used to find all variants relative to the clade-specific reference, and these variants are queried for product names associated with resistance.

        The genes in which there are known resistance-conferring mutations for this pathogen are:

        - Cyp51A
        - HapE
        - COX10 (AFUA_4G08340)
 
        We query `Snippy` results to see if any mutations were identified in those genes. By default, we automatically check for the following loci (which can be overwritten by the user). You will find the mutations next to the locus tag in the `theiaeuk_snippy_variants_hits` column corresponding gene name (see below):

        | **TheiaEuk Search Term** | **Corresponding Gene Name** |
        |---|---|
        | Cyp51A | Cyp51A |
        | HapE | HapE |
        | AFUA_4G08340 | COX10 |

        !!! techdetails "Snippy Variants Technical Details"
            |  | Links |
            | --- | --- |
            | Task | [task_snippy_variants.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snippy_variants.wdl)<br>[task_snippy_gene_query.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snippy_gene_query.wdl) |
            | Software Source Code | [Snippy on GitHub](https://github.com/tseemann/snippy) |
            | Software Documentation | [Snippy on GitHub](https://github.com/tseemann/snippy) |

??? toggle "_Cryptococcus neoformans_"
    When this species is detected by the taxon ID tool an antifungal resistance detection task is deployed.

    ??? task "Snippy Variants: antifungal resistance detection"
        To detect mutations that may confer antifungal resistance, `Snippy` is used to find all variants relative to the clade-specific reference, and these variants are queried for product names associated with resistance.

        The genes in which there are known resistance-conferring mutations for this pathogen are:

        - ERG11 (CNA00300)
        
        We query `Snippy` results to see if any mutations were identified in those genes. By default, we automatically check for the following loci (which can be overwritten by the user). You will find the mutations next to the locus tag in the `theiaeuk_snippy_variants_hits` column corresponding gene name (see below):

        | **TheiaEuk Search Term** | **Corresponding Gene Name** |
        |---|---|
        | CNA00300 | ERG11 |
    
        !!! techdetails "Snippy Variants Technical Details"
            |  | Links |
            | --- | --- |
            | Task | [task_snippy_variants.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snippy_variants.wdl)<br>[task_snippy_gene_query.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snippy_gene_query.wdl) |
            | Software Source Code | [Snippy on GitHub](https://github.com/tseemann/snippy) |
            | Software Documentation | [Snippy on GitHub](https://github.com/tseemann/snippy) |

### Outputs

<div class="searchable-table" markdown="1">

| **Variable** | **Type** | **Description** |
|---|---|---|
| assembly_fasta | File | _De novo_ genome assembly in FASTA format |
| assembly_length | Int | Length of assembly (total number of nucleotides) as determined by QUAST |
| assembler | String | Assembler used in digger_denovo subworkflow |
| assembler_version | String | Version of the assembler used |
| amr_results_csv | File | CSV formatted AMR profile |
| amr_results_pdf | File | PDF formatted AMR profile |
| amr_search_results | File | JSON formatted AMR profile including BLAST results |
| amr_search_docker | String | Docker image used to run AMR_Search |
| amr_search_version | String | Version of AMR_Search libraries used |
| bbduk_docker| String | BBDuk docker image used |
| busco_database | String | BUSCO database used |
| busco_docker | String | BUSCO docker image used |
| busco_report | File | A plain text summary of the results in BUSCO notation |
| busco_results | String | BUSCO results (see above for explanation of BUSCO notation) |
| busco_version | String | BUSCO software version used |
| cg_pipeline_docker | String | Docker file used for running CG-Pipeline on cleaned reads |
| cg_pipeline_report | File | TSV file of read metrics from raw reads, including average read length, number of reads, and estimated genome coverage |
| cladetyper_annotated_reference | String | The annotated reference file for the identified clade, "None" if no clade was identified |
| cladetyper_clade | String | The clade assigned to the input assembly |
| cladetyper_docker_image | String | The Docker container used for the task |
| cladetyper_gambit_version | String | The version of GAMBIT used for the analysis |
| combined_mean_q_clean | Float | Mean quality score for the combined clean reads |
| combined_mean_q_raw | Float | Mean quality score for the combined raw reads |
| combined_mean_readlength_clean | Float | Mean read length for the combined clean reads |
| combined_mean_readlength_raw | Float | Mean read length for the combined raw reads |
| contigs_gfa | File | Assembly graph if spades used for genome assembly |
| est_coverage_clean | Float | Estimated coverage calculated from clean reads and genome length |
| est_coverage_raw | Float | Estimated coverage calculated from raw reads and genome length |
| fastp_html_report | File | The HTML report made with fastp |
| fastp_version | String | Version of fastp software used |
| fastq_scan_clean1_json | File | JSON file output from `fastq-scan` containing summary stats about clean forward read quality and length |
| fastq_scan_clean2_json | File | JSON file output from `fastq-scan` containing summary stats about clean reverse read quality and length |
 fastq_scan_num_reads_clean_pairs | String | Number of read pairs after cleaning as calculated by fastq_scan |
| fastq_scan_num_reads_clean1 | Int | Number of forward reads after cleaning as calculated by fastq_scan |
| fastq_scan_num_reads_clean2 | Int | Number of reverse reads after cleaning as calculated by fastq_scan |
| fastq_scan_num_reads_raw_pairs | String | Number of input read pairs calculated by fastq_scan |
| fastq_scan_num_reads_raw1 | Int | Number of input forward reads calculated by fastq_scan |
| fastq_scan_num_reads_raw2 | Int | Number of input reverse reads calculated by fastq_scan |
| fastq_scan_num_reads_raw_pairs | String | Number of input read pairs calculated by fastq_scan |
| fastq_scan_raw1_json | File | JSON file output from `fastq-scan` containing summary stats about raw forward read quality and length |
| fastq_scan_raw2_json | File | JSON file output from `fastq-scan` containing summary stats about raw reverse read quality and length |
| fastq_scan_version | String | Version of fastq-scan software used |
| fastqc_clean1_html | File | Graphical visualization of clean forward read quality from fastqc to open in an internet browser |
| fastqc_clean2_html | File | Graphical visualization of clean reverse read quality from fastqc to open in an internet browser |
| fastqc_docker | String | Docker container used with fastqc |
| fastqc_num_reads_clean1 | Int | Number of forward reads after cleaning by fastqc |
| fastqc_num_reads_clean2 | Int | Number of reverse reads after cleaning by fastqc |
| fastqc_num_reads_clean_pairs | String | Number of read pairs after cleaning by fastqc |
| fastqc_num_reads_raw1 | Int | Number of input reverse reads by fastqc |
| fastqc_num_reads_raw2 | Int | Number of input reverse reads by fastqc |
| fastqc_num_reads_raw_pairs | String | Number of input read pairs by fastqc |
| fastqc_raw1_html | File | Graphical visualization of raw forward read quality from fastqc to open in an internet browser |
| fastqc_raw2_html | File | Graphical visualization of raw reverse read qualityfrom fastqc to open in an internet browser |
| fastqc_version | String | Version of fastqc software used |
| filtered_contigs_metrics | File | File containing metrics of contigs filtered |
| gambit_closest_genomes | File | CSV file listing genomes in the GAMBIT database that are most similar to the query assembly |
| gambit_db_version | String | Version of GAMBIT used |
| gambit_docker | String | GAMBIT docker file used |
| gambit_predicted_taxon | String | Taxon predicted by GAMBIT |
| gambit_predicted_taxon_rank | String | Taxon rank of GAMBIT taxon prediction |
| gambit_report | File | GAMBIT report in a machine-readable format |
| gambit_version | String | Version of GAMBIT software used |
| n50_value | Int | N50 of assembly calculated by QUAST |
| number_contigs | Int | Total number of contigs in assembly |
| qc_check | String | A string that indicates whether or not the sample passes a set of pre-determined and user-provided QC thresholds |
| qc_standard | File | The user-provided file that contains the QC thresholds used for the QC check |
| quast_gc_percent | Float | The GC percent of your sample |
| quast_report | File | TSV report from QUAST |
| quast_version | String | Software version of QUAST used |
| r1_mean_q_raw | Float | Mean quality score of raw forward reads |
| r1_mean_readlength_raw | Float | Mean read length of raw forward reads |
| r2_mean_q_raw | Float | Mean quality score of raw reverse reads |
| r2_mean_readlength_clean | Float | Mean read length of clean reverse reads |
| rasusa_version | String | Version of rasusa used |
| read1_clean | File | Clean forward reads file |
| read1_subsampled | File | Subsampled read1 file |
| read2_clean | File | Clean reverse reads file |
| read2_subsampled | File | Subsampled read2 file |
| read_screen_raw | String | PASS or FAIL result from raw read screening; FAIL accompanied by the reason(s) for failure |
| read_screen_raw_tsv | File | Raw read screening report TSV depicting read counts, total read base pairs, and estimated genome length |
| read_screen_clean | String | PASS or FAIL result from clean read screening; FAIL accompanied by the reason(s) for failure |
| read_screen_clean_tsv | File | Clean read screening report TSV depicting read counts, total read base pairs, and estimated genome length |
| seq_platform | String | Sequencing platform input by the user |
| theiaeuk_illumina_pe_analysis_date | String | Date of TheiaEuk PE workflow execution |
| theiaeuk_illumina_pe_version | String | TheiaEuk PE workflow version used |
| theiaeuk_snippy_variants_bai | String | BAI file produced by the snippy module |
| theiaeuk_snippy_variants_bam | String | BAM file produced by the snippy module |
| theiaeuk_snippy_variants_coverage_tsv | String | TSV file containing coverage information for each base in the reference genome |
| theiaeuk_snippy_variants_gene_query_results | File | File containing all lines from variants file matching gene query terms |
| theiaeuk_snippy_variants_hits | String | String of all variant file entries matching gene query term |
| theiaeuk_snippy_variants_num_reads_aligned | String | Number of reads aligned by snippy |
| theiaeuk_snippy_variants_num_variants | Int | Number of variants detected by snippy |
| theiaeuk_snippy_variants_outdir_tarball | File | Tar compressed file containing full snippy output directory |
| theiaeuk_snippy_variants_percent_ref_coverage | String | Percent of reference genome covered by snippy |
| theiaeuk_snippy_variants_query | String | The gene query term(s) used to search variant |
| theiaeuk_snippy_variants_query_check | String | Were the gene query terms present in the refence annotated genome file |
| theiaeuk_snippy_variants_reference_genome | File | The reference genome used in the alignment and variant calling |
| theiaeuk_snippy_variants_results | File | The variants file produced by snippy |
| theiaeuk_snippy_variants_summary | File | A file summarizing the variants detected by snippy |
| theiaeuk_snippy_variants_version | String | The version of the snippy_variants module being used |
| trimmomatic_docker | String | Docker image used for trimmomatic |
| trimmomatic_version | String | Version of trimmomatic used |

</div>
