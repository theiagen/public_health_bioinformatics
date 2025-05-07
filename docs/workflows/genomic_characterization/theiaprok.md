# TheiaProk Workflow Series

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Genomic Characterization](../../workflows_overview/workflows_type.md/#genomic-characterization) | [Bacteria](../../workflows_overview/workflows_kingdom.md/#bacteria) | PHB v3.0.1 | Yes, some optional features incompatible | Sample-level |

## TheiaProk Workflows

**The TheiaProk workflows are for the assembly, quality assessment, and characterization of bacterial genomes.** There are currently four TheiaProk workflows designed to accommodate different kinds of input data:

1. Illumina paired-end sequencing **(TheiaProk_Illumina_PE**)
2. Illumina single-end sequencing (**TheiaProk_Illumina_SE)**
3. ONT sequencing (**TheiaProk_ONT**)
4. Genome assemblies (**TheiaProk_FASTA**)

!!! caption "TheiaProk Workflow Diagram"
    ![TheiaProk Workflow Diagram](../../assets/figures/TheiaProk.png)

All input reads are processed through "[core tasks](#core-tasks-performed-for-all-taxa)" in the TheiaProk Illumina and ONT workflows. These undertake read trimming and assembly appropriate to the input data type. TheiaProk workflows subsequently launch default genome characterization modules for quality assessment, species identification, antimicrobial resistance gene detection, sequence typing, and more. **For some taxa identified, "taxa-specific sub-workflows" will be automatically activated, undertaking additional taxa-specific characterization steps.** When setting up each workflow, users may choose to use "optional tasks" as additions or alternatives to tasks run in the workflow by default.

### Inputs

!!! dna ""
    ??? toggle "TheiaProk_Illumina_PE Input Read Data"

        The TheiaProk_Illumina_PE workflow takes in Illumina paired-end read data. Read file names should end with `.fastq` or `.fq`, with the optional addition of `.gz`. When possible, Theiagen recommends zipping files with [gzip](https://www.gnu.org/software/gzip/) before Terra uploads to minimize data upload time.

        By default, the workflow anticipates **2 x 150bp** reads (i.e. the input reads were generated using a 300-cycle sequencing kit). Modifications to the optional parameter for `trim_minlen` may be required to accommodate shorter read data, such as the 2 x 75bp reads generated using a 150-cycle sequencing kit.

    ??? toggle "TheiaProk_Illumina_SE Input Read Data"

        TheiaProk_Illumina_SE takes in Illumina single-end reads. Read file names should end with `.fastq` or `.fq`, with the optional addition of `.gz`. Theiagen highly recommends zipping files with [gzip](https://www.gnu.org/software/gzip/) before uploading to Terra to minimize data upload time & save on storage costs.

        By default, the workflow anticipates **1 x 35 bp** reads  (i.e. the input reads were generated using a 70-cycle sequencing kit). Modifications to the optional parameter for `trim_minlen` may be required to accommodate longer read data.

    ??? toggle "TheiaProk_ONT Input Read Data"

        The TheiaProk_ONT workflow takes in base-called ONT read data. Read file names should end with `.fastq` or `.fq`, with the optional addition of `.gz`. When possible, Theiagen recommends zipping files with [gzip](https://www.gnu.org/software/gzip/) before uploading to Terra to minimize data upload time.

        **The ONT sequencing kit and base-calling approach can produce substantial variability in the amount and quality of read data. Genome assemblies produced by the TheiaProk_ONT workflow must be quality assessed before reporting results.**

    ??? toggle "TheiaProk_FASTA Input Assembly Data"

        The TheiaProk_FASTA workflow takes in assembly files in FASTA format.

<div class="searchable-table" markdown="1">

| **Terra Task name** | **Variable** | **Type** | **Description** | **Default value** | **Terra Status** | **Workflow** |
|---|---|---|---|---|---|---|
| *workflow name | **samplename** | String | Name of sample to be analyzed |  | Required | FASTA, ONT, PE, SE |
| theiaprok_fasta | **assembly_fasta** | File | Assembly file in fasta format |  | Required | FASTA |
| theiaprok_illumina_pe | **read1** | File | Illumina forward read file in FASTQ file format (compression optional) |  | Required | PE |
| theiaprok_illumina_pe | **read2** | File | Illumina reverse read file in FASTQ file format (compression optional) |  | Required | PE |
| theiaprok_illumina_se | **read1** | File | Illumina forward read file in FASTQ file format (compression optional) |  | Required | SE |
| theiaprok_ont | **read1** | File | Base-called ONT read file in FASTQ file format (compression optional) |  | Required | ONT |
| *workflow name | **abricate_db** | String | Database to use with the Abricate tool. Options: NCBI, CARD, ARG-ANNOT, Resfinder, MEGARES, EcOH, PlasmidFinder, Ecoli_VF and VFDB | vfdb | Optional | FASTA, ONT, PE, SE |
| *workflow name | **bakta_db** | String |Database selection for Bakta annotation. Options: "light" (smaller, faster), "full" (more comprehensive), or a Google Storage URI (gs://...) pointing to a custom Bakta database archive (.tar.gz). The selected database will be extracted before annotation. | full | Optional | FASTA, ONT, PE, SE |
| *workflow name | **call_abricate** | Boolean | Set to true to enable the Abricate task | FALSE | Optional | FASTA, ONT, PE, SE |
| *workflow name | **call_ani** | Boolean | Set to true to enable the ANI task | FALSE | Optional | FASTA, ONT, PE, SE |
| *workflow name | **call_kmerfinder** | Boolean | Set to true to enable the kmerfinder task | FALSE | Optional | FASTA, ONT, PE, SE |
| *workflow name | **call_plasmidfinder** | Boolean | Set to true to enable the plasmidfinder task | TRUE | Optional | FASTA, ONT, PE, SE |
| *workflow name | **call_resfinder** | Boolean | Set to true to enable the ResFinder task | FALSE | Optional | FASTA, ONT, PE, SE |
| *workflow name | **city** | String | Will be used in the "city" column in any taxon-specific tables created in the Export Taxon Tables task |  | Optional | FASTA, ONT, PE, SE |
| *workflow name | **collection_date** | String | Will be used in the "collection_date" column in any taxon-specific tables created in the Export Taxon Tables task |  | Optional | FASTA, ONT, PE, SE |
| *workflow name | **county** | String | Will be used in the "county" column in any taxon-specific tables created in the Export Taxon Tables task |  | Optional | FASTA, ONT, PE, SE |
| *workflow name | **expected_taxon** | String | If provided, this input will override the taxonomic assignment made by GAMBIT and launch the relevant taxon-specific submodules. It will also modify the organism flag used by AMRFinderPlus. Example format: "Salmonella enterica" |  | Optional | FASTA, ONT, PE, SE |
| *workflow name | **genome_annotation** | String | If set to "bakta", TheiaProk will use Bakta rather than Prokka to annotate the genome | prokka | Optional | FASTA, ONT, PE, SE |
| *workflow name | **genome_length** | Int | User-specified expected genome length to be used in genome statistics calculations |  | Optional | ONT, PE, SE |
| *workflow name | **max_genome_length** | Int | Maximum genome length able to pass read screening. For TheiaProk_ONT, screening using max_genome_length is skipped by default. | 18040666 | Optional | ONT, PE, SE |
| *workflow name | **min_basepairs** | Int | Minimum number of base pairs able to pass read screening | 2241820 | Optional | ONT, PE, SE |
| *workflow name | **min_coverage** | Int | Minimum genome coverage able to pass read screening. Screening using min_coverage is skipped by default. | 5 | Optional | ONT |
| *workflow name | **min_coverage** | Int | Minimum genome coverage able to pass read screening | 10 | Optional | PE, SE |
| *workflow name | **min_genome_length** | Int | Minimum genome length able to pass read screening. For TheiaProk_ONT, screening using min_genome_length is skipped by default. | 100000 | Optional | ONT, PE, SE |
| *workflow name | **min_proportion** | Int | Minimum proportion of total reads in each read file to pass read screening | 40 | Optional | PE |
| *workflow name | **min_reads** | Int | Minimum number of reads to pass read screening | 5000 | Optional | ONT |
| *workflow name | **min_reads** | Int | Minimum number of reads to pass read screening | 7472 | Optional | PE, SE |
| *workflow name | **originating_lab** | String | Will be used in the "originating_lab" column in any taxon-specific tables created in the Export Taxon Tables task |  | Optional | FASTA, ONT, PE, SE |
| *workflow name | **perform_characterization** | Boolean | Set to "false" if you want to only generate an assembly and relevant QC metrics and skip all characterization tasks | TRUE | Optional | FASTA, ONT, PE, SE |
| *workflow name | **qc_check_table** | File | TSV value with taxons for rows and QC values for columns; internal cells represent user-determined QC thresholds; if provided, turns on the QC Check task.<br>Click on the variable name for an example QC Check table |  | Optional | FASTA, ONT, PE, SE |
| *workflow name | **read1_lane2** | File | If provided, the Concatenate_Illumina_Lanes subworkflow will concatenate all files from the same lane before doing any subsequent analysis |  | Optional | PE, SE |
| *workflow name | **read1_lane3** | File | If provided, the Concatenate_Illumina_Lanes subworkflow will concatenate all files from the same lane before doing any subsequent analysis |  | Optional | PE, SE |
| *workflow name | **read1_lane4** | File | If provided, the Concatenate_Illumina_Lanes subworkflow will concatenate all files from the same lane before doing any subsequent analysis |  | Optional | PE, SE |
| *workflow name | **read2_lane2** | File | If provided, the Concatenate_Illumina_Lanes subworkflow will concatenate all files from the same lane before doing any subsequent analysis |  | Optional | PE, SE |
| *workflow name | **read2_lane3** | File | If provided, the Concatenate_Illumina_Lanes subworkflow will concatenate all files from the same lane before doing any subsequent analysis |  | Optional | PE, SE |
| *workflow name | **read2_lane4** | File | If provided, the Concatenate_Illumina_Lanes subworkflow will concatenate all files from the same lane before doing any subsequent analysis |  | Optional | PE, SE |
| *workflow name | **run_id** | String | Will be used in the "run_id" column in any taxon-specific tables created in the Export Taxon Tables task |  | Optional | FASTA, ONT, PE, SE |
| *workflow name | **seq_method** | String | Will be used in the "seq_id" column in any taxon-specific tables created in the Export Taxon Tables task |  | Optional | FASTA, ONT, PE, SE |
| *workflow name | **skip_mash** | Boolean | If true, skips estimation of genome size and coverage in read screening steps. As a result, providing true also prevents screening using these parameters. | TRUE | Optional | ONT, SE |
| *workflow name | **skip_screen** | Boolean | Option to skip the read screening prior to analysis | FALSE | Optional | ONT, PE, SE |
| *workflow name | **taxon_tables** | File | File indicating data table names to copy samples of a particular taxon to |  | Optional | FASTA, ONT, PE, SE |
| *workflow name | **terra_project** | String | The name of the Terra Project where you want the taxon tables written to |  | Optional | FASTA, ONT, PE, SE |
| *workflow name | **terra_workspace** | String | The name of the Terra Workspace where you want the taxon tables written to |  | Optional | FASTA, ONT, PE, SE |
| *workflow name | **trim_min_length** | Int | Specifies minimum length of each read after trimming to be kept | 25 | Optional | SE |
| *workflow name | **trim_min_length** | Int | Specifies minimum length of each read after trimming to be kept | 75 | Optional | PE |
| *workflow name | **trim_quality_min_score** | Int | Specifies the minimum average quality of bases in a sliding window to be kept | 20 | Optional | PE |
| *workflow name | **trim_quality_trim_score** | Int | Specifies the average quality of bases in a sliding window to be kept | 30 | Optional | SE |
| *workflow name | **trim_window_size** | Int | Specifies window size for trimming (the number of bases to average the quality across) | 4 | Optional | SE |
| *workflow name | **trim_window_size** | Int | Specifies window size for trimming (the number of bases to average the quality across) | 4 | Optional | PE |
| *workflow name | **zip** | String | Will be used in the "zip" column in any taxon-specific tables created in the Export Taxon Tables task |  | Optional | FASTA, ONT, PE, SE |
| abricate | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional | FASTA, ONT, PE, SE |
| abricate | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional | FASTA, ONT, PE, SE |
| abricate | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/abricate:1.0.1-abaum-plasmid | Optional | FASTA, ONT, PE, SE |
| abricate | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional | FASTA, ONT, PE, SE |
| abricate | **min_percent_coverage** | Int | Minimum DNA %coverage for the Abricate task | 80 | Optional | FASTA, ONT, PE, SE |
| abricate | **min_percent_identity** | Int | Minimum DNA %identity for the Abricate task | 80 | Optional | FASTA, ONT, PE, SE |
| amrfinderplus_task | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional | FASTA, ONT, PE, SE |
| amrfinderplus_task | **detailed_drug_class** | Boolean | If set to true, amrfinderplus_amr_classes and amrfinderplus_amr_subclasses outputs will be created | FALSE | Optional | FASTA, ONT, PE, SE |
| amrfinderplus_task | **disk_size** | Boolean | Amount of storage (in GB) to allocate to the AMRFinderPlus task | 50 | Optional | FASTA, ONT, PE, SE |
| amrfinderplus_task | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/ncbi-amrfinderplus:4.0.19-2024-12-18.1 | Optional | FASTA, ONT, PE, SE |
| amrfinderplus_task | **hide_point_mutations** | Boolean | If set to true, point mutations are not reported | FALSE | Optional | FASTA, ONT, PE, SE |
| amrfinderplus_task | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional | FASTA, ONT, PE, SE |
| amrfinderplus_task | **min_percent_coverage** | Float | Minimum proportion of reference gene covered for a BLAST-based hit (Methods BLAST or PARTIAL)." Attribute should be a float ranging from 0-1, such as 0.6 (equal to 60% coverage) | 0.5 | Optional| FASTA, ONT, PE, SE |
| amrfinderplus_task | **min_percent_identity** | Float | "Minimum identity for a blast-based hit hit (Methods BLAST or PARTIAL). -1 means use a curated threshold if it exists and 0.9 otherwise. Setting this value to something other than -1 will override any curated similarity cutoffs." Attribute should be a float ranging from 0-1, such as 0.95 (equal to 95% identity) | 0.9 | Optional | FASTA, ONT, PE, SE |
| amrfinderplus_task | **separate_betalactam_genes**  | Boolean | Report beta-Lactam AMR genes separated out by all beta-lactam and the respective beta-lactam subclasses | FALSE | Optional | FASTA, ONT, PE, SE |
| ani | **ani_threshold** | Float | ANI value threshold must be surpassed in order to output the ani_top_species_match. If a genome does not surpass this threshold (and the percent_bases_aligned_threshold) then the ani_top_species_match output String will show a warning instead of a genus & species. | 80 | Optional | FASTA, ONT, PE, SE |
| ani | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional | FASTA, ONT, PE, SE |
| ani | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional | FASTA, ONT, PE, SE |
| ani | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/staphb/mummer:4.0.0-rgdv2 | Optional | FASTA, ONT, PE, SE |
| ani | **mash_filter** | Float | Mash distance threshold over which ANI is not calculated  | 0.9 | Optional | FASTA, ONT, PE, SE |
| ani | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional | FASTA, ONT, PE, SE |
| ani | **percent_bases_aligned_threshold** | Float | Threshold regarding the proportion of bases aligned between the query genome and reference genome. If a genome does not surpass this threshold (and the ani_threshold) then the ani_top_species_match output String will show a warning instead of a genus & species. | 70 | Optional | FASTA, ONT, PE, SE |
| ani | **ref_genome** | File | If not set, uses all 43 genomes in RGDv2 |  | Optional | FASTA, ONT, PE, SE |
| bakta | **bakta_opts** | String | Parameters to pass to bakta from <https://github.com/oschwengers/bakta#usage> |  | Optional | FASTA, ONT, PE, SE |
| bakta | **compliant** | Boolean | If true, forces Genbank/ENA/DDJB compliance | FALSE | Optional | FASTA, ONT, PE, SE |
| bakta | **cpu** | Int | Number of CPUs to allocate to the task | 8 | Optional | FASTA, ONT, PE, SE |
| bakta | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional | FASTA, ONT, PE, SE |
| bakta | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/biocontainers/bakta:1.5.1--pyhdfd78af_0 | Optional | FASTA, ONT, PE, SE |
| bakta | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 16 | Optional | FASTA, ONT, PE, SE |
| bakta | **prodigal_tf** | File | Prodigal training file to use for CDS prediction by bakta |  | Optional | FASTA, ONT, PE, SE |
| bakta | **proteins** | Boolean |  | FALSE | Optional | FASTA, ONT, PE, SE |
| busco | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional | FASTA, ONT, PE, SE |
| busco | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional | FASTA, ONT, PE, SE |
| busco | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/ezlabgva/busco:v5.7.1_cv1 | Optional | FASTA, ONT, PE, SE |
| busco | **eukaryote** | Boolean | Assesses eukaryotic organisms, rather than prokaryotic organisms | FALSE | Optional | FASTA, ONT, PE, SE |
| busco | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional | FASTA, ONT, PE, SE |
| cg_pipeline_clean | **cg_pipe_opts** | String | Options to pass to CG-Pipeline for clean read assessment  | --fast | Optional | PE, SE |
| cg_pipeline_clean | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional | PE, SE |
| cg_pipeline_clean | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional | PE, SE |
| cg_pipeline_clean | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/lyveset:1.1.4f | Optional | PE, SE |
| cg_pipeline_clean | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional | PE, SE |
| cg_pipeline_clean | **read2** | File | Internal component, do not modify |  | Do not modify, Optional | SE |
| cg_pipeline_raw | **cg_pipe_opts** | String | Options to pass to CG-Pipeline for raw read assessment  | --fast | Optional | PE, SE |
| cg_pipeline_raw | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional | PE, SE |
| cg_pipeline_raw | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional | PE, SE |
| cg_pipeline_raw | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/lyveset:1.1.4f | Optional | PE, SE |
| cg_pipeline_raw | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional | PE, SE |
| cg_pipeline_raw | **read2** | File | Internal component, do not modify |  | Do not modify, Optional | SE |
| clean_check_reads | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional | ONT, PE, SE |
| clean_check_reads | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional | ONT, PE, SE |
| clean_check_reads | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/bactopia/gather_samples:2.0.2 | Optional | ONT, PE, SE |
| clean_check_reads | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional | ONT, PE, SE |
| clean_check_reads | **organism** | String | Internal component, do not modify |  | Do not modify, Optional | ONT, PE, SE |
| clean_check_reads | **workflow_series** | String | Internal component, do not modify |  | Do not modify, Optional | ONT, PE, SE |
| digger_denovo | **assembler** | String | Assembler to use (spades, skesa, megahit) | spades | Optional | PE, SE |
| digger_denovo | **assembly_options** | String | String | Assembler-specific options that you might choose for the selected assembler | | Optional | PE, SE |
| digger_denovo | **bwa_cpu** | Int | Number of CPU cores for BWA alignment | 6 | Optional | PE, SE |
| digger_denovo | **bwa_disk_size** | Int | Disk space in GB for BWA alignment | 100 | Optional | PE, SE |
| digger_denovo | **bwa_docker** | String | Docker image for BWA alignment | us-docker.pkg.dev/general-theiagen/staphb/ivar:1.3.1-titan | Optional | PE, SE |
| digger_denovo | **bwa_memory** | Int | Memory in GB for BWA alignment | 16 | Optional | PE, SE |
| digger_denovo | **filter_contigs_cpu** | Int | Number of CPU cores for contig filtering | 1 | Optional | PE, SE |
| digger_denovo | **filter_contigs_disk_size** | Int | Disk space in GB for contig filtering | 50 | Optional | PE, SE |
| digger_denovo | **filter_contigs_docker** | String | Docker image for contig filtering | us-docker.pkg.dev/general-theiagen/theiagen/shovilter:0.2 | Optional | PE, SE |
| digger_denovo | **filter_contigs_memory** | Int | Memory in GB for contig filtering | 8 | Optional | PE, SE |
| digger_denovo | **filter_contigs_min_coverage** | Float | Minimum coverage threshold for contig filtering | 2.0 | Optional | PE, SE |
| digger_denovo | **filter_contigs_min_length** | Int | Minimum length threshold for contig filtering | 200 | Optional | PE, SE |
| digger_denovo | **filter_contigs_skip_coverage_filter** | Boolean | Skip filtering contigs based on coverage | false | Optional | PE, SE |
| digger_denovo | **filter_contigs_skip_homopolymer_filter** | Boolean | Skip filtering contigs containing homopolymers | false | Optional | PE, SE |
| digger_denovo | **filter_contigs_skip_length_filter** | Boolean | Skip filtering contigs based on length | false | Optional | PE, SE |
| digger_denovo | **kmers** | String | K-mer sizes for assembly (comma-separated) | | Optional | PE, SE |
| digger_denovo | **megahit_cpu** | Int | Number of CPU cores for MEGAHIT assembler | 4 | Optional | PE, SE |
| digger_denovo | **megahit_disk_size** | Int | Disk space in GB for MEGAHIT assembler | 100 | Optional | PE, SE |
| digger_denovo | **megahit_docker** | String | Docker image for MEGAHIT assembler | us-docker.pkg.dev/general-theiagen/theiagen/megahit:1.2.9 | Optional | PE, SE |
| digger_denovo | **megahit_memory** | Int | Memory in GB for MEGAHIT assembler | 16 | Optional | PE, SE |
| digger_denovo | **min_contig_length** | Int | Minimum contig length to retain in final assembly | 200 | Optional | PE, SE |
| digger_denovo | **pilon_cpu** | Int | Number of CPU cores for Pilon polishing | 8 | Optional | PE, SE |
| digger_denovo | **pilon_disk_size** | Int | Disk space in GB for Pilon polishing | 100 | Optional | PE, SE |
| digger_denovo | **pilon_docker** | String | Docker image for Pilon polishing | us-docker.pkg.dev/general-theiagen/biocontainers/pilon:1.24--hdfd78af_0 | Optional | PE, SE |
| digger_denovo | **pilon_memory** | Int | Memory in GB for Pilon polishing | 32 | Optional | PE, SE |
| digger_denovo | **run_filter_contigs** | Boolean | Whether to run contig filtering step | true | Optional | PE, SE |
| digger_denovo | **skesa_cpu** | Int | Number of CPU cores for SKESA assembler | 4 | Optional | PE, SE |
| digger_denovo | **skesa_disk_size** | Int | Disk space in GB for SKESA assembler | 50 | Optional | PE, SE |
| digger_denovo | **skesa_docker** | String | Docker image for SKESA assembler | us-docker.pkg.dev/general-theiagen/staphb/skesa:2.4.0 | Optional | PE, SE |
| digger_denovo | **skesa_memory** | Int | Memory in GB for SKESA assembler | 4 | Optional | PE, SE |
| digger_denovo | **spades_cpu** | Int | Number of CPU cores for SPAdes assembler | 16 | Optional | PE, SE |
| digger_denovo | **spades_disk_size** | Int | Disk space in GB for SPAdes assembler | 100 | Optional | PE, SE |
| digger_denovo | **spades_docker** | String | Docker image for SPAdes assembler | us-docker.pkg.dev/general-theiagen/staphb/spades:4.1.0 | Optional | PE, SE |
| digger_denovo | **spades_memory** | Int | Memory in GB for SPAdes assembler | 32 | Optional | PE, SE |
| digger_denovo | **spades_type** | String | SPAdes assembly mode (isolate, meta, rna, etc.), more can be found [here](https://ablab.github.io/spades/running.html) | isolate | Optional | PE, SE |
| digger_denovo | **use_pilon** | Boolean | Whether to run Pilon polishing after assembly | false | Optional | PE, SE |
| flye_denovo | **auto_medaka_model** | Boolean | If true, medaka will automatically select the best Medaka model for assembly | TRUE | Optional | ONT |
| flye_denovo | **bandage_cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional | ONT |
| flye_denovo | **bandage_disk_size** | Int | Amount of storage (in GB) to allocate to the task | 10 | Optional | ONT |
| flye_denovo | **bandage_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 4 | Optional | ONT |
| flye_denovo | **dnaapler_cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional | ONT |
| flye_denovo | **dnaapler_disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional | ONT |
| flye_denovo | **dnaapler_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 16 | Optional | ONT |
| flye_denovo | **dnaapler_mode** | String | Dnaapler-specific inputs | all | Optional | ONT |
| flye_denovo | **filtercontigs_cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional | ONT |
| flye_denovo | **filtercontigs_disk_size** | Int | Amount of storage (in GB) to allocate to the task | 10 | Optional | ONT |
| flye_denovo | **filtercontigs_min_length** | Int | Minimum contig length to keep | 1000 | Optional | ONT |
| flye_denovo | **filtercontigs_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 16 | Optional | ONT |
| flye_denovo    | **flye_read_type**   | String | Specifies the type of sequencing reads. Options: `--nano-raw` (default), `--nano-corr`, `--nano-hq`, `--pacbio-raw`, `--pacbio-corr`, `--pacbio-hifi`. Refer to Flye documentation for details on each type. | `--nano-raw`        | Optional | `--nano-raw`   |
| flye_denovo | **flye_genome_length** | Int | User-specified expected genome length to be used in genome statistics calculations |  | Optional | ONT |
| flye_denovo | **flye_asm_coverage** | Int | Reduced coverage for initial disjointig assembly |  | Optional | ONT |
| flye_denovo | **flye_polishing_iterations** | Int | Default polishing iterations | 1 | Optional | ONT |
| flye_denovo | **flye_minimum_overlap** | Int | Minimum overlap between reads |  | Optional | ONT |
| flye_denovo | **flye_read_error_rate** | Float | Maximum expected read error rate |  | Optional | ONT |
| flye_denovo | **flye_uneven_coverage_mode** | Boolean |  | FALSE | Optional | ONT |
| flye_denovo | **flye_keep_haplotypes** | Boolean | If true keep haplotypes | FALSE | Optional | ONT |
| flye_denovo | **flye_no_alt_contigs** | Boolean | If true, do not generate alternative contigs | FALSE | Optional | ONT |
| flye_denovo | **flye_scaffold** | Boolean | If true, scaffold | FALSE | Optional | ONT |
| flye_denovo | **flye_additional_parameters** | String | Any extra Flye-specific parameters |  | Optional | ONT |
| flye_denovo | **flye_cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional | ONT |
| flye_denovo | **flye_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 32 | Optional | ONT |
| flye_denovo | **flye_disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional | ONT |
| flye_denovo | **illumina_read1** | File | If Illumina reads are provided, flye_denovo subworkflow will perform Illumina polishing |  | Optional | ONT |
| flye_denovo | **illumina_read2** | File | If Illumina reads are provided, flye_denovo subworflow will perform Illumina polishing |  | Optional | ONT |
| flye_denovo | **medaka_model** | String | In order to obtain the best results, the appropriate model must be set to match the sequencer's basecaller model; this string takes the format of {pore}_{device}_{caller variant}_{caller_version}. See also <https://github.com/nanoporetech/medaka?tab=readme-ov-file#models>. If this is being run on legacy data it is likely to be r941_min_hac_g507.  |   r1041_e82_400bps_sup_v5.0.0 | Optional | ONT |
| flye_denovo | **medaka_cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional | ONT |
| flye_denovo | **medaka_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 16 | Optional | ONT |
| flye_denovo | **medaka_disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional | ONT |
| flye_denovo | **polisher** | String | The polishing tool to use for assembly | medaka | Optional | ONT |
| flye_denovo | **polishing_rounds** | Int | The number of polishing rounds to conduct for medaka or racon (without Illumina) | 1 | Optional | ONT |
| flye_denovo | **porechop_cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional | ONT |
| flye_denovo | **porechop_disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional | ONT |
| flye_denovo | **porechop_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 16 | Optional | ONT |
| flye_denovo | **porechop_trimopts** | String | Options to pass to Porechop for trimming |  | Optional | ONT |
| flye_denovo | **polypolish_pair_orientation** | String | Polypolish-specific inputs |  | Optional | ONT |
| flye_denovo | **polypolish_low_percentile_threshold** | Float | Polypolish-specific inputs |  | Optional | ONT |
| flye_denovo | **polypolish_high_percentile_threshold** | Float | Polypolish-specific inputs |  | Optional | ONT |
| flye_denovo | **polypolish_fraction_invalid** | Float | Polypolish-specific inputs |  | Optional | ONT |
| flye_denovo | **polypolish_fraction_valid** | Float | Polypolish-specific inputs |  | Optional | ONT |
| flye_denovo | **polypolish_maximum_errors** | Int | Polypolish-specific inputs |  | Optional | ONT |
| flye_denovo | **polypolish_minimum_depth** | Int | Polypolish-specific inputs |  | Optional | ONT |
| flye_denovo | **polypolish_careful** | Boolean | Polypolish-specific inputs | FALSE | Optional | ONT |
| flye_denovo | **polypolish_cpu** | Int | Polypolish cpu | 1 | Optional | ONT |
| flye_denovo | **polypolish_memory** | Int | Polypolish memory | 8 | Optional | ONT |
| flye_denovo | **polypolish_disk_size** | Int | Polypolish disk size | 100 | Optional | ONT |
| flye_denovo | **racon_cpu** | Int | Number of CPUs to allocate to the task | 8 | Optional | ONT |
| flye_denovo | **racon_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 16 | Optional | ONT |
| flye_denovo | **racon_disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional | ONT |
| flye_denovo | **read1** | File | ONT read file in FASTQ file format (compression optional) |  | Optional | ONT |
| flye_denovo | **run_porechop** | Boolean | If true, trims reads before assembly using Porechop | FALSE | Optional | ONT |
| flye_denovo | **skip_polishing** | Boolean | If true, skips polishing | FALSE | Optional | ONT |
| export_taxon_tables | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional | FASTA, ONT, PE, SE |
| export_taxon_tables | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 25 | Optional | FASTA, ONT, PE, SE |
| export_taxon_tables | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-06-21 | Optional | FASTA, ONT, PE, SE |
| export_taxon_tables | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional | FASTA, ONT, PE, SE |
| gambit | **cpu** | Int | Number of CPUs to allocate to the task | 8 | Optional | FASTA, ONT, PE, SE |
| gambit | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional | FASTA, ONT, PE, SE |
| gambit | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/gambit:1.0.0 | Optional | FASTA, ONT, PE, SE |
| gambit | **gambit_db_genomes** | File | User-provided database of assembled query genomes; requires complementary signatures file. If not provided, uses default database, "/gambit-db" | gs://gambit-databases-rp/2.0.0/gambit-metadata-2.0.0-20240628.gdb | Optional | FASTA, ONT, PE, SE |
| gambit | **gambit_db_signatures** | File | User-provided signatures file; requires complementary genomes file. If not specified, the file from the docker container will be used.  | gs://gambit-databases-rp/2.0.0/gambit-signatures-2.0.0-20240628.gs | Optional | FASTA, ONT, PE, SE |
| gambit | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 16 | Optional | FASTA, ONT, PE, SE |
| kmerfinder | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional | FASTA, ONT, PE, SE |
| kmerfinder | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional | FASTA, ONT, PE, SE |
| kmerfinder | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/biocontainers/kmerfinder:3.0.2--hdfd78af_0 | Optional | FASTA, ONT, PE, SE |
| kmerfinder | **kmerfinder_args** | String | Kmerfinder additional arguments |  | Optional | FASTA, ONT, PE, SE |
| kmerfinder | **kmerfinder_db** | String | Bacterial database for KmerFinder | gs://theiagen-public-files-rp/terra/theiaprok-files/kmerfinder_bacteria_20230911.tar.gz | Optional | FASTA, ONT, PE, SE |
| kmerfinder | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 32 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **abricate_abaum_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/abricate:1.0.1-abaum-plasmid | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **abricate_abaum_min_percent_coverage** | Int | Minimum DNA percent coverage |  | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **abricate_abaum_min_percent_identity** | Int | Minimum DNA percent identity; set to 95 because there is a strict threshold of 95% identity for typing purposes | 95 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **abricate_vibrio_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/abricate:1.0.1-abaum-plasmid | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **abricate_vibrio_min_percent_coverage** | Int | Minimum DNA percent coverage  | 80 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **abricate_vibrio_min_percent_identity** | Int | Minimum DNA percent identity | 80 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **agrvate_agr_typing_only** | Boolean | Set to true to skip agr operon extraction and frameshift detection | False  | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **agrvate_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/biocontainers/agrvate:1.0.2--hdfd78af_0 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **amr_search_cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **amr_search_disk_size** | Int | Amount of storage (in GB) to allocate to the task | 50 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **amr_search_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/amrsearch:0.2.1 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **amr_search_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **assembly_only** | Boolean | Internal component, do not modify |  | Do not modify, Optional | ONT, PE, SE |
| merlin_magic | **call_poppunk** | Boolean | If "true", runs PopPUNK for GPSC cluster designation for S. pneumoniae | TRUE | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **call_shigeifinder_reads_input** | Boolean | If set to "true", the ShigEiFinder task will run again but using read files as input instead of the assembly file. Input is shown but not used for TheiaProk_FASTA. | FALSE | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **call_stxtyper** | Boolean | If set to "true", the StxTyper task will run on all samples regardless of the `gambit_predicted_taxon` output. Useful if you suspect a non-E.coli or non-Shigella sample contains stx genes. | FALSE | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **call_tbp_parser** | Boolean | If set to "true", activates the tbp_parser module and results in more outputs, including tbp_parser_looker_report_csv, tbp_parser_laboratorian_report_csv,  tbp_parser_lims_report_csv, tbp_parser_coverage_report, and tbp_parser_genome_percent_coverage | FALSE | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **cauris_cladetyper_docker_image** | String | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **cladetyper_kmer_size** | Int | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **cladetyper_ref_clade1** | File | *Provide an empty file if running TheiaProk on the command-line |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **cladetyper_ref_clade1_annotated** | File | *Provide an empty file if running TheiaProk on the command-line |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **cladetyper_ref_clade2** | File | *Provide an empty file if running TheiaProk on the command-line |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **cladetyper_ref_clade2_annotated** | File | *Provide an empty file if running TheiaProk on the command-line |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **cladetyper_ref_clade3** | File | *Provide an empty file if running TheiaProk on the command-line |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **cladetyper_ref_clade3_annotated** | File | *Provide an empty file if running TheiaProk on the command-line |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **cladetyper_ref_clade4** | File | *Provide an empty file if running TheiaProk on the command-line |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **cladetyper_ref_clade4_annotated** | File | *Provide an empty file if running TheiaProk on the command-line |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **cladetyper_ref_clade5** | File | *Provide an empty file if running TheiaProk on the command-line |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **cladetyper_ref_clade5_annotated** | File | *Provide an empty file if running TheiaProk on the command-line |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **clockwork_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/cdcgov/varpipe_wgs_with_refs:2bc7234074bd53d9e92a1048b0485763cd9bbf6f4d12d5a1cc82bfec8ca7d75e | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **ectyper_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/biocontainers/ectyper:1.0.0--pyhdfd78af_1 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **ectyper_h_min_percent_coverage** | Int | Minumum percent coverage required for an H antigen allele match | 50 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **ectyper_h_min_percent_identity** | Int | Percent identity required for an H antigen allele match | 95 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **ectyper_o_min_percent_coverage** | Int | Minumum percent coverage required for an O antigen allele match | 90 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **ectyper_o_min_percent_identity** | Int | Percent identity required for an O antigen allele match | 90 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **ectyper_print_alleles** | Boolean | Set to true to print the allele sequences as the final column | False  | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **ectyper_verify** | Boolean | Set to true to enable E. coli species verification | False  | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **emmtypingtool_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/emmtypingtool:0.0.1 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **genotyphi_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/mykrobe:0.11.0 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **hicap_min_broken_gene_percent_identity** | Float | Minimum percentage identity to consider a broken gene | 0.80 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **hicap_broken_gene_length** | Int | Minimum length to consider a broken gene | 60 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **hicap_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/biocontainers/hicap:1.0.3--py_0 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **hicap_min_gene_percent_coverage** | Float | Minimum percentage coverage to consider a single gene complete | 0.80 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **hicap_min_gene_percent_identity** | Float | Minimum percentage identity to consider a single gene complete | 0.70 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **kaptive_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/kaptive:2.0.3 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **kaptive_low_gene_percent_identity** | Float | Percent identity threshold for what counts as a low identity match in the gene BLAST search | 95 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **kaptive_min_coverage** | Float | Minimum required percent identity for the gene BLAST search via tBLASTn | 80 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **kaptive_min_percent_identity** | Float | Minimum required percent coverage for the gene BLAST search via tBLASTn | 90 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **kaptive_start_end_margin** | Int | Determines flexibility in identifying the start and end of a locus - if this value is 10, a locus match that is missing the first 8 base pairs will still count as capturing the start of the locus | 10 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **kleborate_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/kleborate:2.2.0 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **kleborate_min_percent_coverage** | Float | Minimum alignment percent coverage for main results | 80 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **kleborate_min_percent_identity** | Float | Minimum alignment percent identity for main results | 90 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **kleborate_min_kaptive_confidence** | String | {None,Low,Good,High,Very_high,Perfect} Minimum Kaptive confidence to call K/O loci - confidence levels below this will be reported as unknown | Good | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **kleborate_min_spurious_percent_coverage** | Float | Minimum alignment percent coverage for spurious results | 40 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **kleborate_min_spurious_percent_identity** | Float | Minimum alignment percent identity for spurious results | 80 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **kleborate_skip_kaptive** | Boolean | Equivalent to --kaptive_k --kaptive_ | False  | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **kleborate_skip_resistance** | Boolean | Set to true to turn on resistance genes screening (default: no resistance gene screening) | False  | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **legsta_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/biocontainers/legsta:0.5.1--hdfd78af_2 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **lissero_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/biocontainers/lissero:0.4.9--py_0 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **lissero_min_percent_coverage** | Float | Minimum percent coverage of the gene to accept a match | 95 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **lissero_min_percent_identity** | Float | Minimum percent identity to accept a match | 95 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **meningotype_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/biocontainers/meningotype:0.8.5--pyhdfd78af_0 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **ngmaster_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/ngmaster:1.0.0 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **ont_data** | Boolean | Internal component, do not modify |  | Do not modify, Optional | FASTA, PE, SE |
| merlin_magic | **paired_end** | Boolean | Internal component, do not modify |  | Do not modify, Optional | ONT, PE |
| merlin_magic | **pasty_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/pasty:1.0.3 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **pasty_min_percent_coverage** | Int | Minimum coverage of a O-antigen to be considered for serogrouping by pasty | 95 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **pasty_min_percent_identity** | Int | Minimum percent identity for a blast hit to be considered for serogrouping | 95 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **pbptyper_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/pbptyper:1.0.4 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **pbptyper_min_percent_coverage** | Int | Minimum percent coverage to count a hit | 90 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **pbptyper_min_percent_identity** | Int | Minimum percent identity to count a hit | 90 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **poppunk_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/poppunk:2.4.0 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **poppunk_gps_clusters_csv** | File | Poppunk database file *Provide an empty or local file if running TheiaProk on the command-line | gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6_clusters.csv | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **poppunk_gps_dists_npy** | File | Poppunk database file *Provide an empty or local file if running TheiaProk on the command-line | gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6.dists.npy | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **poppunk_gps_dists_pkl** | File | Poppunk database file *Provide an empty or local file if running TheiaProk on the command-line | gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6.dists.pkl | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **poppunk_gps_external_clusters_csv** | File | Poppunk database file *Provide an empty or local file if running TheiaProk on the command-line | gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6_external_clusters.csv | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **poppunk_gps_fit_npz** | File | Poppunk database file *Provide an empty or local file if running TheiaProk on the command-line | gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6_fit.npz | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **poppunk_gps_fit_pkl** | File | Poppunk database file *Provide an empty or local file if running TheiaProk on the command-line | gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6_fit.pkl | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **poppunk_gps_graph_gt** | File | Poppunk database file *Provide an empty or local file if running TheiaProk on the command-line | gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6_graph.gt | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **poppunk_gps_h5** | File | Poppunk database file *Provide an empty or local file if running TheiaProk on the command-line | gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6.h5 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **poppunk_gps_qcreport_txt** | File | Poppunk database file *Provide an empty or local file if running TheiaProk on the command-line | gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6_qcreport.txt | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **poppunk_gps_refs** | File | Poppunk database file *Provide an empty or local file if running TheiaProk on the command-line | gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6.refs | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **poppunk_gps_refs_dists_npy** | File | Poppunk database file *Provide an empty or local file if running TheiaProk on the command-line | gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6.refs.dists.npy | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **poppunk_gps_refs_dists_pkl** | File | Poppunk database file *Provide an empty or local file if running TheiaProk on the command-line | gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6.refs.dists.pkl | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **poppunk_gps_refs_graph_gt** | File | Poppunk database file *Provide an empty or local file if running TheiaProk on the command-line | gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6refs_graph.gt | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **poppunk_gps_refs_h5** | File | Poppunk database file *Provide an empty or local file if running TheiaProk on the command-line | gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6.refs.h5 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **poppunk_gps_unword_clusters_csv** | File | Poppunk database file *Provide an empty or local file if running TheiaProk on the command-line | gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6_unword_clusters.csv | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **read1** | File | Internal component, do not modify |  | Do not modify, Optional | FASTA |
| merlin_magic | **read2** | File | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, SE |
| merlin_magic | **run_amr_search** | Boolean | If set to true AMR_Search workflow will be run if species is part of supported taxon, see AMR_Search docs. | False | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **seqsero2_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/seqsero2:1.2.1 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **seroba_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/seroba:1.0.2 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **serotypefinder_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/serotypefinder:2.0.1 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **shigatyper_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/shigatyper:2.0.5 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **shigeifinder_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/shigeifinder:1.3.5 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **sistr_cpu** | Int | The number of CPU cores to allocate for the task | 8 | Optional | FASTA, ONT, PE, SE  |
| merlin_magic | **sistr_disk_size** | Int | The disk size (in GB) to allocate for the task  | 100 | Optional | FASTA, ONT, PE, SE  |
| merlin_magic | **sistr_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/biocontainers/sistr_cmd:1.1.1--pyh864c0ab_2 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **sistr_memory** | Int | The amount of memory (in GB) to allocate for the task. | 32 | Optional | FASTA, ONT, PE, SE  |
| merlin_magic | **sistr_use_full_cgmlst_db** | Boolean | Set to true to use the full set of cgMLST alleles which can include highly similar alleles. By default the smaller "centroid" alleles or representative alleles are used for each marker | False  | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **snippy_base_quality** | Int | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **snippy_gene_query_docker_image** | String | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **snippy_map_qual** | Int | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **snippy_maxsoft** | Int | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **snippy_min_coverage** | Int | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **snippy_min_frac** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **snippy_min_quality** | Int | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **snippy_query_gene** | String | Internal component, do not modify |  | Do not modify, Optional | FASTA, PE, SE |
| merlin_magic | **snippy_reference_afumigatus** | File | *Provide an empty file if running TheiaProk on the command-line |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **snippy_reference_calbicans** | File | *Provide an empty file if running TheiaProk on the command-line |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **snippy_reference_cryptoneo** | File | *Provide an empty file if running TheiaProk on the command-line |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **snippy_variants_docker_image** | String | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **sonneityping_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/mykrobe:0.12.1 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **sonneityping_mykrobe_opts** | String | Additional options for mykrobe in sonneityping |  | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **spatyper_do_enrich** | Boolean | Set to true to enable PCR product enrichment | False  | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **spatyper_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/biocontainers/spatyper:0.3.3--pyhdfd78af_3 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **srst2_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/srst2:0.2.0-vcholerae | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **srst2_gene_max_mismatch** | Int | Maximum number of mismatches for SRST2 to call a gene as present | 2000 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **srst2_max_divergence** | Int | Maximum divergence, in percentage, for SRST2 to call a gene as present | 20 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **srst2_min_coverage** | Int | Minimum breadth of coverage for SRST2 to call a gene as present | 80 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **srst2_min_depth** | Int | Minimum depth of coverage for SRST2 to call a gene as present  | 5 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **srst2_min_edge_depth** | Int | Minimum edge depth for SRST2 to call a gene as present  | 2 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **stxtyper_cpu** | Int | The number of CPU cores to allocate for the task. | 1 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **stxtyper_disk_size** | Int | Amount of storage (in GB) to allocate to the task | 50 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **stxtyper_docker_image** | String | The Docker container to use for the task | `us-docker.pkg.dev/general-theiagen/staphb/stxtyper:1.0.24` | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **stxtyper_enable_debug** | Boolean | When enabled, additional messages are printed and files in `$TMPDIR` are not removed after running | FALSE | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **stxtyper_memory** | Int | Amount of memory (in GB) to allocate to the task | 4 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **staphopia_sccmec_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/biocontainers/staphopia-sccmec:1.0.0--hdfd78af_0 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbp_parser_config** | File | The configuration file to use, in YAML format (overrides all other arguments except input_json and input_bam) |  | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbp_parser_add_cs_lims** | Boolean | Set to true add cycloserine results to the LIMS report | FALSE | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbp_parser_coverage_regions_bed** | File | A bed file that lists the regions to be considered for QC |  | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbp_parser_min_percent_coverage** | Int | The minimum coverage for a region to pass QC in tbp_parser | 100 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbp_parser_debug** | Boolean | Activate the debug mode on tbp_parser; increases logging outputs | TRUE | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbp_parser_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/tbp-parser:2.4.5 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbp_parser_etha237_frequency** | Float | Minimum frequency for a mutation in ethA at protein position 237 to pass QC in tbp-parser | 0.1 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbp_parser_expert_rule_regions_bed** | File | A file that contains the regions where R mutations and expert rules are applied |  | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbp_parser_min_depth** | Int | Minimum depth for a variant to pass QC in tbp_parser | 10 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbp_parser_min_frequency** | Int | The minimum frequency for a mutation to pass QC | 0.1 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbp_parser_min_read_support** | Int | The minimum read support for a mutation to pass QC | 10 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbp_parser_operator** | String | Fills the "operator" field in the tbp_parser output files | Operator not provided | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbp_parser_output_seq_method_type** | String | Fills out the "seq_method" field in the tbp_parser output files | Sequencing method not provided | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbp_parser_rpob449_frequency** | Float | Minimum frequency for a mutation at protein position 449 to pass QC in tbp-parser | 0.1 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbp_parser_rrl_frequency** | Float | Minimum frequency for a mutation in rrl to pass QC in tbp-parser | 0.1 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbp_parser_rrl_read_support** | Int | Minimum read support for a mutation in rrl to pass QC in tbp-parser | 10 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbp_parser_rrs_frequency** | Float | Minimum frequency for a mutation in rrs to pass QC in tbp-parser | 0.1 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbp_parser_rrs_read_support** | Int | Minimum read support for a mutation in rrs to pass QC in tbp-parser | 10 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbp_parser_tngs_data** | Boolean | Set to true to enable tNGS-specific parameters and runs in tbp-parser | FALSE | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbprofiler_custom_db** | File | TBProfiler uses by default the TBDB database; if you have a custom database you wish to use, you must provide a custom database in this field and set tbprofiler_run_custom_db to true |  | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbprofiler_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/tbprofiler:6.6.3 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbprofiler_mapper** | String | The mapping tool used in TBProfiler to align the reads to the reference genome; see TBProfiler’s original documentation for available options. | bwa | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbprofiler_min_af** | Float | The minimum allele frequency to call a variant | 0.1 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbprofiler_min_depth** | Int | The minimum depth for a variant to be called. | 10 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbprofiler_run_cdph_db** | Boolean | TBProfiler uses by default the TBDB database; set this value to "true" to use the WHO v2 database with customizations for CDPH | FALSE | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbprofiler_run_custom_db** | Boolean | TBProfiler uses by default the TBDB database; if you have a custom database you wish to use, you must set this value to true and provide a custom database in the tbprofiler_custom_db field | FALSE | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbprofiler_variant_caller** | String | Select a different variant caller for TBProfiler to use by writing it in this block; see TBProfiler’s original documentation for available options. | GATK | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **tbprofiler_variant_calling_params** | String | Enter additional variant calling parameters in this free text input to customize how the variant caller works in TBProfiler | | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **theiaeuk** | Boolean | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| merlin_magic | **vibecheck_lineage_barcodes** | String   | Feather formatted lineage barcodes to use instead of default O1 barcodes |  | Optional | PE |
| merlin_magic | **vibecheck_subsampling_fraction** | Float | Fraction of reads to use in classification. | 0.2  | Optional | PE |
| merlin_magic | **vibecheck_skip_subsampling**  | Boolean | When enabled, will not subsample reads prior to classification. Will increase computation time | False | Optional | PE |
| merlin_magic | **vibecheck_docker_image** | String | The Docker container to use for the task | watronfire/vibecheck:2025.02.24 | Optional | PE |
| merlin_magic | **virulencefinder_min_percent_coverage** | Float | The threshold for minimum coverage |  | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **virulencefinder_database** | String | The specific database to use | virulence_ecoli | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **virulencefinder_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/virulencefinder:2.0.4 | Optional | FASTA, ONT, PE, SE |
| merlin_magic | **virulencefinder_min_percent_identity** | Float | The threshold for minimum blast identity |  | Optional | FASTA, ONT, PE, SE |
| nanoplot_clean | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional | ONT |
| nanoplot_clean | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional | ONT |
| nanoplot_clean | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/nanoplot:1.40.0 | Optional | ONT |
| nanoplot_clean | **max_length** | Int | Maximum read length for nanoplot | 100000 | Optional | ONT |
| nanoplot_clean | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 16 | Optional | ONT |
| nanoplot_raw | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional | ONT |
| nanoplot_raw | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional | ONT |
| nanoplot_raw | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/nanoplot:1.40.0 | Optional | ONT |
| nanoplot_raw | **max_length** | Int | Maximum read length for nanoplot | 100000 | Optional | ONT |
| nanoplot_raw | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 16 | Optional | ONT |
| plasmidfinder | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional | FASTA, ONT, PE, SE |
| plasmidfinder | **database** | String | User-specified database |  | Optional | FASTA, ONT, PE, SE |
| plasmidfinder | **database_path** | String | Path to user-specified database |  | Optional | FASTA, ONT, PE, SE |
| plasmidfinder | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 50 | Optional | FASTA, ONT, PE, SE |
| plasmidfinder | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/plasmidfinder:2.1.6 | Optional | FASTA, ONT, PE, SE |
| plasmidfinder | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional | FASTA, ONT, PE, SE |
| plasmidfinder | **method_path** | String | Path to files for a user-specified method to use (blast or kma) |  | Optional | FASTA, ONT, PE, SE |
| plasmidfinder | **min_percent_coverage** | Float | Threshold for minimum coverage, default threshold from PlasmidFinder CLI tool is used (0.60) | 0.6 | Optional | FASTA, ONT, PE, SE |
| plasmidfinder | **min_percent_identity** | Float | Threshold for mininum blast identity, default threshold from PlasmidFinder CLI tool is used (0.90). This default differs from the default of the PlasmidFinder webtool (0.95) | 0.9 | Optional | FASTA, ONT, PE, SE |
| prokka | **compliant** | Boolean | Forces Genbank/ENA/DDJB compliant headers in Prokka output files | TRUE | Optional | FASTA, ONT, PE, SE |
| prokka | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional | FASTA, ONT, PE, SE |
| prokka | **disk_size** | String | Amount of storage (in GB) to allocate to the PlasmidFinder task | 100 | Optional | FASTA, ONT, PE, SE |
| prokka | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/prokka:1.14.5 | Optional | FASTA, ONT, PE, SE |
| prokka | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 16 | Optional | FASTA, ONT, PE, SE |
| prokka | **prodigal_tf** | File | <https://github.com/tseemann/prokka#option---prodigaltf> |  | Optional | FASTA, ONT, PE, SE |
| prokka | **prokka_arguments** | String | Any additional <https://github.com/tseemann/prokka#command-line-options> |  | Optional | FASTA, ONT, PE, SE |
| prokka | **proteins** | Boolean | FASTA file of trusted proteins for Prokka to first use for annotations | FALSE | Optional | FASTA, ONT, PE, SE |
| qc_check_task | **assembly_length_unambiguous** | Int | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| qc_check_task | **assembly_mean_coverage** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| qc_check_task | **combined_mean_q_clean** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, SE |
| qc_check_task | **combined_mean_q_raw** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, SE |
| qc_check_task | **combined_mean_readlength_clean** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, SE |
| qc_check_task | **combined_mean_readlength_raw** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, SE |
| qc_check_task | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional | FASTA, ONT, PE, SE |
| qc_check_task | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional | FASTA, ONT, PE, SE |
| qc_check_task | **docker** | String | The Docker container to use for the task |  "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16" | Optional | FASTA, ONT, PE, SE |
| qc_check_task | **est_coverage_clean** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA |
| qc_check_task | **est_coverage_raw** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA |
| qc_check_task | **kraken_human** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| qc_check_task | **kraken_human_dehosted** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| qc_check_task | **kraken_sc2** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| qc_check_task | **kraken_sc2_dehosted** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| qc_check_task | **kraken_target_organism** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| qc_check_task | **kraken_target_organism_dehosted** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| qc_check_task | **meanbaseq_trim** | String | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| qc_check_task | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional | FASTA, ONT, PE, SE |
| qc_check_task | **midas_secondary_genus_abundance** | Int | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT |
| qc_check_task | **midas_secondary_genus_coverage** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT |
| qc_check_task | **num_reads_clean1** | Int | Internal component, do not modify |  | Do not modify, Optional | FASTA |
| qc_check_task | **num_reads_clean2** | Int | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, SE |
| qc_check_task | **num_reads_raw1** | Int | Internal component, do not modify |  | Do not modify, Optional | FASTA |
| qc_check_task | **num_reads_raw2** | Int | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, SE |
| qc_check_task | **number_Degenerate** | Int | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| qc_check_task | **number_N** | Int | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| qc_check_task | **percent_reference_coverage** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| qc_check_task | **r1_mean_q_clean** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA |
| qc_check_task | **r1_mean_q_raw** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA |
| qc_check_task | **r1_mean_readlength_clean** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA |
| qc_check_task | **r1_mean_readlength_raw** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA |
| qc_check_task | **r2_mean_q_clean** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, SE |
| qc_check_task | **r2_mean_q_raw** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, SE |
| qc_check_task | **r2_mean_readlength_clean** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, SE |
| qc_check_task | **r2_mean_readlength_raw** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, SE |
| qc_check_task | **sc2_s_gene_mean_coverage** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| qc_check_task | **sc2_s_gene_percent_coverage** | Float | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| qc_check_task | **vadr_num_alerts** | String | Internal component, do not modify |  | Do not modify, Optional | FASTA, ONT, PE, SE |
| quast | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional | FASTA, ONT, PE, SE |
| quast | **disk_size** | String | Amount of storage (in GB) to allocate to the Quast task | 100 | Optional | FASTA, ONT, PE, SE |
| quast | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/quast:5.0.2 | Optional | FASTA, ONT, PE, SE |
| quast | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional | FASTA, ONT, PE, SE |
| quast | **min_contig_length** | Int | Lower threshold for a contig length in bp. Shorter contigs won’t be taken into account | 500 | Optional | FASTA, ONT, PE, SE |
| raw_check_reads | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional | ONT, PE, SE |
| raw_check_reads | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional | ONT, PE, SE |
| raw_check_reads | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/bactopia/gather_samples:2.0.2 | Optional | ONT, PE, SE |
| raw_check_reads | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional | ONT, PE, SE |
| raw_check_reads | **organism** | String | Internal component, do not modify |  | Do not modify, Optional | ONT, PE, SE |
| raw_check_reads | **workflow_series** | String | Internal component, do not modify |  | Do not modify, Optional | ONT, PE, SE |
| read_QC_trim | **adapters** | File | A file containing the sequence of the adapters used during library preparation, used in the BBDuk task |  | Optional | PE, SE |
| read_QC_trim | **artic_guppyplex_cpu** | Int | Internal component, do not modify| 8 | Optional | ONT |
| read_QC_trim | **artic_guppyplex_disk_size** | Int | Internal component, do not modify| 100 | Optional | ONT |
| read_QC_trim | **artic_guppyplex_docker** | String | Internal component, do not modify| us-docker.pkg.dev/general-theiagen/staphb/artic-ncov2019:1.3.0-medaka-1.4.3 | Optional | ONT |
| read_QC_trim | **artic_guppyplex_memory** | Int | Internal component, do not modify| 16 | Optional | ONT |
| read_QC_trim | **bbduk_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional | PE, SE |
| read_QC_trim | **call_kraken** | Boolean | Set to true to launch Kraken2; if true, you must provide a kraken_db | FALSE | Optional | ONT, PE, SE |
| read_QC_trim | **call_midas** | Boolean | Set to true to launch Midas | TRUE | Optional | PE, SE |
| read_QC_trim | **downsampling_coverage** | Float | The depth to downsample to with Rasusa | 150 | Optional | ONT |
| read_QC_trim | **fastp_args** | String | Additional arguments to pass to fastp | -g -5 20 -3 20 | Optional | SE |
| read_QC_trim | **fastp_args** | String | Additional arguments to pass to fastp | "--detect_adapter_for_pe -g -5 20 -3 20 | Optional | PE |
| read_QC_trim | **kraken_cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional | ONT, PE, SE |
| read_QC_trim | **kraken_db** | File | Kraken2 database file; must be provided in call_kraken is true |  | Optional | ONT, PE, SE |
| read_QC_trim | **kraken_disk_size** | Int | GB of storage to request for VM used to run the kraken2 task. Increase this when using large (>30GB kraken2 databases such as the "k2_standard" database) | 100 | Optional | ONT, PE, SE |
| read_QC_trim | **kraken_docker_image** | Int | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/staphb/kraken2:2.1.2-no-db" | Optional | ONT |
| read_QC_trim | **kraken_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional | ONT, PE, SE |
| read_QC_trim | **max_length** | Int | Internal component, do not modify |  | Do not modify, Optional | ONT |
| read_QC_trim | **min_length** | Int | Internal component, do not modify |  | Do not modify, Optional | ONT |
| read_QC_trim | **midas_db** | File | Midas database file | gs://theiagen-large-public-files-rp/terra/theiaprok-files/midas/midas_db_v1.2.tar.gz | Optional | PE, SE |
| read_QC_trim | **phix** | File | A file containing the phix used during Illumina sequencing; used in the BBDuk task |  | Optional | PE, SE |
| read_QC_trim | **nanoq_cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional | ONT |
| read_QC_trim | **nanoq_disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional | ONT |
| read_QC_trim | **nanoq_docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/biocontainers/nanoq:0.9.0--hec16e2b_1" | Optional | ONT |
| read_QC_trim | **nanoq_max_read_length** | Int | The maximum read length to keep after trimming | 100000 | Optional | ONT |
| read_QC_trim | **nanoq_max_read_qual** | Int | The maximum read quality to keep after trimming | 40 | Optional | ONT |
| read_QC_trim | **nanoq_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional | ONT |
| read_QC_trim | **nanoq_min_read_length** | Int | The minimum read length to keep after trimming | 500 | Optional | ONT |
| read_QC_trim | **nanoq_min_read_qual** | Int | The minimum read quality to keep after trimming | 10 | Optional | ONT |
| read_QC_trim | **ncbi_scrub_cpu** | Int | Internal component, do not modify| 4 | Optional | ONT |
| read_QC_trim | **ncbi_scrub_disk_size** | Int | Internal component, do not modify| 100 | Optional | ONT |
| read_QC_trim | **ncbi_scrub_docker** | String | Internal component, do not modify| "us-docker.pkg.dev/general-theiagen/ncbi/sra-human-scrubber:2.2.1" | Optional | ONT |
| read_QC_trim | **ncbi_scrub_memory** | Int | Internal component, do not modify| 8 | Optional | ONT |
| read_QC_trim | **rasusa_bases** | String | Explicitly set the number of bases required e.g., 4.3kb, 7Tb, 9000, 4.1MB. If this option is given, --coverage and --genome-size are ignored | | Optional | ONT |
| read_QC_trim | **rasusa_cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional | ONT |
| read_QC_trim | **rasusa_disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional | ONT |
| read_QC_trim | **rasusa_docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/staphb/rasusa:2.1.0" | Optional | ONT |
| read_QC_trim | **rasusa_fraction_of_reads** | Float | Subsample to a fraction of the reads - e.g., 0.5 samples half the reads | | Optional | ONT |
| read_QC_trim | **rasusa_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional | ONT |
| read_QC_trim | **rasusa_number_of_reads** | Int | Subsample to a specific number of reads | | Optional | ONT |
| read_QC_trim | **rasusa_seed** | Int | Random seed to use | | Optional | ONT |
| read_QC_trim | **read_processing** | String | Read trimming software to use, either "trimmomatic" or "fastp" | trimmomatic | Optional | PE, SE |
| read_QC_trim | **read_qc** | String | Allows the user to decide between fastq_scan (default) and fastqc for the evaluation of read quality. | fastq_scan | Optional | PE, SE |
| read_QC_trim | **run_prefix** | String | Internal component, do not modify |  | Do not modify, Optional | ONT |
| read_QC_trim | **target_organism** | String | This string is searched for in the kraken2 outputs to extract the read percentage |  | Optional | ONT, PE, SE |
| read_QC_trim | **trimmomatic_args** | String | Additional arguments to pass to trimmomatic. "-phred33" specifies the Phred Q score encoding which is almost always phred33 with modern sequence data. | -phred33 | Optional | PE, SE |
| resfinder_task | **acquired** | Boolean | Set to true to tell ResFinder to identify acquired resistance genes | TRUE | Optional | FASTA, ONT, PE, SE |
| resfinder_task | **call_pointfinder** | Boolean | Set to true to enable detection of point mutations. | FALSE | Optional | FASTA, ONT, PE, SE |
| resfinder_task | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional | FASTA, ONT, PE, SE |
| resfinder_task | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional | FASTA, ONT, PE, SE |
| resfinder_task | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/resfinder:4.1.11 | Optional | FASTA, ONT, PE, SE |
| resfinder_task | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional | FASTA, ONT, PE, SE |
| resfinder_task | **min_percent_coverage** | Float | Minimum coverage breadth of a gene for it to be identified | 0.5 | Optional | FASTA, ONT, PE, SE |
| resfinder_task | **min_percent_identity** | Float | Minimum identity for ResFinder to identify a gene | 0.9 | Optional | FASTA, ONT, PE, SE |
| ts_mlst | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional | FASTA, ONT, PE, SE |
| ts_mlst | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 50 | Optional | FASTA, ONT, PE, SE |
| ts_mlst | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/mlst:2.23.0-2024-12-31 | Optional | FASTA, ONT, PE, SE |
| ts_mlst | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional | FASTA, ONT, PE, SE |
| ts_mlst | **min_percent_coverage** | Float | Minimum % breadth of coverage to report an MLST allele | 10 | Optional | FASTA, ONT, PE, SE |
| ts_mlst | **min_percent_identity** | Float | Minimum % identity to known MLST gene to report an MLST allele | 95 | Optional | FASTA, ONT, PE, SE |
| ts_mlst | **minscore** | Float | Minimum <https://github.com/tseemann/mlst#scoring-system> to assign an MLST profile | 50 | Optional | FASTA, ONT, PE, SE |
| ts_mlst | **nopath** | Boolean | true = use mlst --nopath. If set to false, filename paths are not stripped from FILE column in output TSV | TRUE | Optional | FASTA, ONT, PE, SE |
| ts_mlst | **scheme** | String | Don’t autodetect the MLST scheme; force this scheme on all inputs (see <https://github.com/tseemann/mlst/blob/master/db/scheme_species_map.tab> for accepted strings) | None | Optional | FASTA, ONT, PE, SE |
| version_capture | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional | FASTA, ONT, PE, SE |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) | | Optional | FASTA, ONT, PE, SE |
</div>

!!! tip "Skip Characterization"
    Ever wanted to skip characterization? Now you can! Set the optional input `perform_characterization` to **`false`** to only generate an assembly and run assembly QC.

### Core Tasks (performed for all taxa)

??? task "`versioning`: Version Capture for TheiaProk"

    The `versioning` task captures the workflow version from the GitHub (code repository) version.
        
    !!! techdetails "Version Capture Technical details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_versioning.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/task_versioning.wdl) |

??? task "`concatenate_illumina_lanes`: Concatenate Multi-Lane Illumina FASTQs ==_for Illumina only_=="

    The `concatenate_illumina_lanes` task concatenates Illumina FASTQ files from multiple lanes into a single file. This task only runs if the `read1_lane2` input file has been provided. All read1 lanes are concatenated together and are used in subsequent tasks, as are the read2 lanes. These concatenated files are also provided as output.

    !!! techdetails "Concatenate Illumina Lanes Technical Details"
        The `concatenate_illumina_lanes` task is run before any downstream steps take place.
        
        |  | Links |
        | --- | --- |
        | Task | [wf_concatenate_illumina_lanes.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/file_handling/wf_concatenate_illumina_lanes.wdl)

??? task "`screen`: Total Raw Read Quantification and Genome Size Estimation"

    The [`screen`](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/comparisons/task_screen.wdl) task ensures the quantity of sequence data is sufficient to undertake genomic analysis. It uses [`fastq-scan`](https://github.com/rpetit3/fastq-scan) and bash commands for quantification of reads and base pairs, and [mash](https://mash.readthedocs.io/en/latest/index.html) sketching to estimate the genome size and its coverage. At each step, the results are assessed relative to pass/fail criteria and thresholds that may be defined by optional user inputs. Samples are run through all threshold checks, regardless of failures, and the workflow will terminate after the `screen` task if any thresholds are not met:

    1. Total number of reads: A sample will fail the read screening task if its total number of reads is less than or equal to `min_reads`.
    2. The proportion of basepairs reads in the forward and reverse read files: A sample will fail the read screening if fewer than `min_proportion` basepairs are in either the reads1 or read2 files.
    3. Number of basepairs: A sample will fail the read screening if there are fewer than `min_basepairs` basepairs
    4. Estimated genome size:  A sample will fail the read screening if the estimated genome size is smaller than `min_genome_size` or bigger than `max_genome_size`.
    5. Estimated genome coverage: A sample will fail the read screening if the estimated genome coverage is less than the `min_coverage`.

    Read screening is undertaken on both the raw and cleaned reads. The task may be skipped by setting the `skip_screen` variable to true.

    Default values vary between the PE and SE workflow. The rationale for these default values can be found below. If two default values are shown, the first is for Illumina workflows and the second is for ONT.

    | Variable  | Default Value | Rationale |
    | --- | --- | --- |
    | `skip_screen` | false | Set to false to avoid waste of compute resources processing insufficient data |
    | `min_reads` | 7472 or 5000 | Calculated from the minimum number of base pairs required for 20x coverage of Nasuia deltocephalinicola genome, the smallest known bacterial genome as of 2019-08-07 (112,091 bp), divided by 300 (the longest Illumina read length) or 5000 (estimate of ONT read length) |
    | `min_basepairs` | 2241820 | Should be greater than 20x coverage of Nasuia deltocephalinicola, the smallest known bacterial genome (112,091 bp) |
    | `min_genome_length` | 100000 | Based on the Nasuia deltocephalinicola genome - the smallest known bacterial genome (112,091 bp) |
    | `max_genome_length` | 18040666 | Based on the Minicystis rosea genome, the biggest known bacterial genome (16,040,666 bp), plus an additional 2 Mbp to cater for potential extra genomic material |
    | `min_coverage` | 10 or 5 | A bare-minimum average per base coverage across the genome required for genome characterization. Note, a higher per base coverage coverage would be required for high-quality phylogenetics. |
    | `min_proportion` | 40 | Neither read1 nor read2 files should have less than 40% of the total number of reads. For paired-end data only |
    
    !!! techdetails "Screen Technical Details"    
        There is a single WDL task for read screening that contains two separate sub-tasks, one used for PE data and the other for SE data. The `screen` task is run twice, once for raw reads and once for clean reads.
        
        |  | TheiaProk_Illumina_PE | TheiaProk_Illumina_SE and TheiaProk_ONT |
        | --- | --- | --- |
        | Task | [task_screen.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/comparisons/task_screen.wdl#L3) (PE sub-task) | [task_screen.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/comparisons/task_screen.wdl#L147) (SE sub-task) |

#### Illumina Data Core Tasks

??? task "`read_QC_trim`: Read Quality Trimming, Adapter Removal, Quantification, and Identification"

    `read_QC_trim` is a sub-workflow within TheiaMeta that removes low-quality reads, low-quality regions of reads, and sequencing adapters to improve data quality. It uses a number of tasks, described below.

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

    **Key Components of the MIDAS Database**

    1. **Species Groups**: 
    - MIDAS clusters bacterial genomes based on 96.5% sequence identity, forming over 5,950 species groups from 31,007 genomes. These groups align with the gold-standard species definition (95% ANI), ensuring highly accurate species identification.

    2. **Genomic Data Structure**:
    - **Marker Genes**: Contains 15 universal single-copy genes used to estimate species abundance.
    - **Representative Genome**: Each species group has a selected representative genome, which minimizes genetic variation and aids in accurate SNP identification.
    - **Pan-genome**: The database includes clusters of non-redundant genes, with options for multi-level clustering (e.g., 99%, 95%, 90% identity), enabling MIDAS to identify gene content within strains at various clustering thresholds.

    3. **Taxonomic Annotation**: 
    - Genomes are annotated based on consensus Latin names. Discrepancies in name assignments may occur due to factors like unclassified genomes or genus-level ambiguities.

    ---

    **Using the Default MIDAS Database**

    TheiaProk uses the pre-loaded MIDAS database in Terra (see input table for current version) by default for bacterial species detection in metagenomic data, requiring no additional setup.

    **How to Set Up the Default MIDAS Database**

    Users can also build their own custom MIDAS database if they want to include specific genomes or configurations. This custom database can replace the default MIDAS database used in Terra. To build a custom MIDAS database, follow the [MIDAS GitHub guide on building a custom database](https://github.com/snayfach/MIDAS/blob/master/docs/build_db.md). Once the database is built, users can upload it to a Google Cloud Storage bucket or Terra workkspace and provide the link to the database in the `midas_db` input variable.

    Alternatively to `MIDAS`, the `Kraken2` task can also be turned on through setting the `call_kraken` input variable as `true` for the identification of reads to detect contamination with non-target taxa.

    Kraken2 is a bioinformatics tool originally designed for metagenomic applications. It has additionally proven valuable for validating taxonomic assignments and checking contamination of single-species (e.g. bacterial isolate) whole genome sequence data. A database must be provided if this optional module is activated, through the kraken_db optional input. A list of suggested databases can be found on [Kraken2 standalone documentation](../standalone/kraken2.md).

    !!! techdetails "read_QC_trim Technical Details"
                
        |  | Links |
        | --- | --- |
        | Sub-workflow | [wf_read_QC_trim_pe.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_read_QC_trim_pe.wdl)<br>[wf_read_QC_trim_se.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_read_QC_trim_se.wdl) |
        | Tasks | [task_fastp.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_fastp.wdl)<br>[task_trimmomatic.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_trimmomatic.wdl)<br>[task_bbduk.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_bbduk.wdl)<br>[task_fastq_scan.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_fastq_scan.wdl)<br>[task_midas.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/contamination/task_midas.wdl)<br>[task_kraken2.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/contamination/task_kraken2.wdl)|
        | Software Source Code | [fastp](https://github.com/OpenGene/fastp); [Trimmomatic](https://github.com/usadellab/Trimmomatic); [fastq-scan](https://github.com/rpetit3/fastq-scan); [MIDAS](https://github.com/snayfach/MIDAS); [Kraken2](https://github.com/DerrickWood/kraken2)|
        | Software Documentation | [fastp](https://github.com/OpenGene/fastp); [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic); [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/); [fastq-scan](https://github.com/rpetit3/fastq-scan); [MIDAS](https://github.com/snayfach/MIDAS); [Kraken2](https://github.com/DerrickWood/kraken2/wiki) |
        | Original Publication(s) | [Trimmomatic: a flexible trimmer for Illumina sequence data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4103590/)<br>[fastp: an ultra-fast all-in-one FASTQ preprocessor](https://academic.oup.com/bioinformatics/article/34/17/i884/5093234?login=false)<br>[An integrated metagenomics pipeline for strain profiling reveals novel patterns of bacterial transmission and biogeography](https://pubmed.ncbi.nlm.nih.gov/27803195/)<br>[Improved metagenomic analysis with Kraken 2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0) |

??? task "`CG-Pipeline`: Assessment of Read Quality, and Estimation of Genome Coverage"

    The`cg_pipeline` task generates metrics about read quality and estimates the coverage of the genome using the "run_assembly_readMetrics.pl" script from [CG-Pipeline](https://github.com/lskatz/CG-Pipeline/). The genome coverage estimates are calculated using both using raw and cleaned reads, using either a user-provided `genome_size` or the estimated genome length generated by QUAST.

    !!! techdetails "CG-Pipeline Technical Details"
        The `cg_pipeline` task is run twice in TheiaProk, once with raw reads, and once with clean reads.
        
        |  | Links |
        | --- | --- |
        | Task | [task_cg_pipeline.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_cg_pipeline.wdl) |
        | Software Source Code | [CG-Pipeline on GitHub](https://github.com/lskatz/CG-Pipeline/) |
        | Software Documentation | [CG-Pipeline on GitHub](https://github.com/lskatz/CG-Pipeline/) |
        | Original Publication(s) | [A computational genomics pipeline for prokaryotic sequencing projects](https://academic.oup.com/bioinformatics/article/26/15/1819/188418) |

??? task "`digger_denovo`: _De novo_ Assembly"

    De Novo assembly will be undertaken only for samples that have sufficient read quantity and quality, as determined by the `screen` task assessment of clean reads. 

    In TheiaProk, assembly is performed using the [digger_denovo](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_digger_denovo.wdl), which is a hat tip to [Shovill](https://github.com/tseemann/shovill) pipeline. This undertakes the assembly with one of three assemblers ([SKESA](https://github.com/ncbi/SKESA) (default), [SPAdes](https://github.com/ablab/spades), [Megahit](https://github.com/voutcn/megahit)), but also performs a number of post processing steps for assembly polishing and contig filtering. Pilon can optionally be run if `use_pilon` is set to true. On defualt, the contig filtering task is set to run, which will remove any homopolymers, contigs below a specificied length, and contigs with coverage below a specified minimum coverage. This can be turned off by setting `run_filter_contigs` to `false`. 
    
    ??? toggle "What is _de novo_  assembly?"
        _De novo_  assembly is the process or product of attempting to reconstruct a genome from scratch (without prior knowledge of the genome) using sequence reads. Assembly of fungal genomes from short-reads will produce multiple contigs per chromosome rather than a single contiguous sequence for each chromosome.
        
    !!! techdetails "Digger-Denovo Technical Details"
        |  | Links |
        | --- | --- |
        | TheiaEuk WDL SubWorkflow | [wf_digger_denovo.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_digger_denovo.wd) |
        | Software Source Code | [digger_denovo on GitHub](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_digger_denovo.wd) |

#### ONT Data Core Tasks

??? task "`read_QC_trim_ont`: Read Quality Trimming, Quantification, and Identification"

    `read_QC_trim_ont` is a sub-workflow within TheiaProk_ONT that filters low-quality reads and trims low-quality regions of reads. It uses several tasks, described below.

    **Estimated genome length**:

    By default, an estimated genome length is set to 5 Mb, which is around 0.7 Mb higher than the average bacterial genome length, according to the information collated [here](https://github.com/CDCgov/phoenix/blob/717d19c19338373fc0f89eba30757fe5cfb3e18a/assets/databases/NCBI_Assembly_stats_20240124.txt). This estimate can be overwritten by the user, and is used by `RASUSA`.

    **Plotting and quantifying long-read sequencing data:** `nanoplot`

    Nanoplot is used for the determination of mean quality scores, read lengths, and number of reads. This task is run once with raw reads as input and once with clean reads as input. If QC has been performed correctly, you should expect **fewer** clean reads than raw reads.

    **Read subsampling:** Samples are automatically randomly subsampled to 150X coverage using `RASUSA`.

    **Plasmid prediction:** Plasmids are identified using replicon sequences used for typing from [PlasmidFinder](https://cge.food.dtu.dk/services/PlasmidFinder/).

    **Read filtering:** Reads are filtered by length and quality using `nanoq`. By default, sequences with less than 500 basepairs and quality score lower than 10 are filtered out to improve assembly accuracy.

    !!! techdetails "read_QC_trim_ont Technical Details"
        
        TheiaProk_ONT calls a sub-workflow listed below, which then calls the individual tasks:
        
        | Workflow | **TheiaProk_ONT** |
        | --- | --- |
        | Sub-workflow | [wf_read_QC_trim_ont.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_read_QC_trim_ont.wdl) |
        | Tasks | [task_nanoplot.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_nanoplot.wdl) [task_fastq_scan.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_fastq_scan.wdl) [task_rasusa.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/task_rasusa.wdl) [task_nanoq.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_nanoq.wdl)
        | Software Source Code | [fastq-scan](https://github.com/rpetit3/fastq-scan), [NanoPlot](https://github.com/wdecoster/NanoPlot), [RASUSA](https://github.com/mbhall88/rasusa), [nanoq](https://github.com/esteinig/nanoq) |
        | Original Publication(s) | [NanoPlot paper](https://academic.oup.com/bioinformatics/article/39/5/btad311/7160911)<br>[RASUSA paper](https://doi.org/10.21105/joss.03941)<br>[Nanoq Paper](https://doi.org/10.21105/joss.02991)<br> |

??? task "`Flye`: _De novo_ Assembly"

    `flye_denovo` is a sub-workflow that performs _de novo_ assembly using Flye for ONT data and supports additional polishing and visualization steps.
    
    !!! tip "Ensure correct medaka model is selected if performing medaka polishing"
        In order to obtain the best results, the appropriate model must be set to match the sequencer's basecaller model; this string takes the format of {pore}\_{device}\_{caller variant}\_{caller_version}. See also <https://github.com/nanoporetech/medaka?tab=readme-ov-file#models>. If `flye` is being run on legacy data the medaka model will likely be `r941_min_hac_g507`. Recently generated data will likely be suited by the default model of `r1041_e82_400bps_sup_v5.0.0`.

    The detailed steps and tasks are as follows:

    ??? toggle "`Porechop`: Read Trimming (optional; off by default)"
        Read trimming is optional and can be enabled by setting the `run_porchop` input variable to true.

        Porechop is a tool for finding and removing adapters from ONT data. Adapters on the ends of reads are trimmed, and when a read has an adapter in the middle, the read is split into two.

        !!! techdetails "Porechop Technical Details"
            |  | Links |
            | --- | --- |
            | WDL Task | [task_porechop.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_porechop.wdl) |
            | Software Source Code | [Porechop on GitHub](https://github.com/rrwick/Porechop) |
            | Software Documentation | [https://github.com/rrwick/Porechop#porechop](https://github.com/rrwick/Porechop#porechop) |

    ??? toggle "`Flye`: _De novo_ Assembly"
        Flye is a _de novo_ assembler for long read data using repeat graphs. Compared to de Bruijn graphs, which require exact k-mer matches, repeat graphs can use approximate matches which better tolerates the error rate of ONT data.

        `flye_read_type` specifies the type of sequencing reads being used for assembly. This parameter significantly impacts the assembly process and should match the characteristics of your input data. Below are the available options:

        - `--nano-hq` (default): Optimized for ONT high-quality reads, such as Guppy5+ SUP or Q20 (<5% error). Recommended for ONT reads processed with Guppy5 or newer
        - `--nano-raw`: For ONT regular reads, pre-Guppy5 (<20% error)
        - `--nano-corr`: ONT reads corrected with other methods (<3% error)
        - `--pacbio-raw`: PacBio regular CLR reads (<20% error)
        - `--pacbio-corr`: PacBio reads corrected with other methods (<3% error)
        - `--pacbio-hifi`: PacBio HiFi reads (<1% error)
        
        Refer to the Flye documentation for detailed guidance on selecting the appropriate `flye_read_type` based on your sequencing data and additional optional paramaters.

        !!! techdetails "Flye Technical Details"
            |  | Links |
            | --- | --- |
            | WDL Task | [task_flye.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_flye.wdl) |
            | Software Source Code | [Flye on GitHub](https://github.com/fenderglass/Flye) |
            | Software Documentation | [Flye Documentation](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md) |
            | Original Publication(s) | [Assembly of long, error-prone reads using repeat graphs](https://www.nature.com/articles/s41587-019-0072-8) |

    ??? toggle "`Bandage`: Graph Visualization"
        Bandage creates _de novo_ assembly graphs containing the assembled contigs and the connections between those contigs.

        !!! techdetails "Bandage Technical Details"
            |  | Links |
            | --- | --- |
            | WDL Task | [task_bandage_plot.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_bandage_plot.wdl) |
            | Software Source Code | [Bandage on GitHub](https://github.com/rrwick/Bandage) |
            | Software Documentation | [Bandage Documentation](https://github.com/rrwick/Bandage#bandage) |
            | Original Publication(s) | [Bandage: interactive visualization of _de novo_ genome assemblies](https://academic.oup.com/bioinformatics/article/31/20/3350/196114) |

    ??? toggle "`Polypolish`: Hybrid Assembly Polishing ==_for ONT and Illumina data_=="
        If short reads are provided with the optional `illumina_read1` and `illumina_read2` inputs, Polypolish will use those short-reads to correct errors in the long-read assemblies. Uniquely, Polypolish uses the short-read alignments where each read is aligned to _all_ possible locations, meaning that repeat regions will have error correction.
    
        !!! techdetails "Polypolish Technical Details"
            |  | Links |
            | --- | --- |
            | Task | [task_polypolish.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/polishing/task_polypolish.wdl) |
            | Software Source Code | [Polypolish on GitHub](https://github.com/rrwick/Polypolish) |
            | Software Documentation | [Polypolish Documentation](https://github.com/rrwick/Polypolish#polypolish) |
            | Original Publication(s) | [Polypolish: short-read polishing of long-read bacterial genome assemblies](https://doi.org/10.1371/journal.pcbi.1009802)<br>[How low can you go? Short-read polishing of Oxford Nanopore bacterial genome assemblies](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001254) |

    ??? toggle "`Medaka`: Polishing of Flye assembly (default; optional)"
        Polishing is optional and can be skipped by setting the `skip_polishing` variable to true. If polishing is skipped, then neither Medaka or Racon will run.

        Medaka is the default assembly polisher used in TheiaProk. Racon may be used alternatively, and if so, Medaka will not run. Medaka uses the raw reads to polish the assembly and generate a consensus sequence. 

        Importantly, Medaka requires knowing the model that was used to generate the read data. There are several ways to provide this information:

        - Automatic Model Selection: Automatically determines the most appropriate Medaka model based on the input data, ensuring optimal polishing results without manual intervention. 
        - User-Specified Model Override: Allows users to specify a particular `Medaka model` if automatic selection does not yield the desired outcome or for specialized use cases.
        - Default Model: If both automatic model selection fails and no user-specified model is provided, Medaka defaults to the predefined fallback model `r1041_e82_400bps_sup_v5.0.0`. 

        !!! info "Medaka Model Resolution Process" 
            Medaka's automatic model selection uses the `medaka tools resolve_model` command to identify the appropriate model for polishing. This process relies on metadata embedded in the input file, which is typically generated by the basecaller. If the automatic selection fails to identify a suitable model, Medaka gracefully falls back to the default model to maintain workflow continuity. **Users should verify the chosen model and consider specifying a model override if necessary.**

        !!! techdetails "Medaka Technical Details"
            |  | Links |
            | --- | --- |
            | WDL Task | [task_medaka.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/polishing/task_medaka.wdl) |
            | Software Source Code | [Medaka on GitHub](https://github.com/nanoporetech/medaka) |
            | Software Documentation | [Medaka Documentation](https://github.com/nanoporetech/medaka#medaka) |

    ??? toggle "`Racon`: Polishing of Flye assembly (alternative; optional)"
        Polishing is optional and can be skipped by setting the `skip_polishing` variable to true. If polishing is skipped, then neither Medaka or Racon will run.

        `Racon` is an alternative to using `medaka` for assembly polishing, and can be run by setting the `polisher` input to "racon".  Racon is a consensus algorithm designed for refining raw de novo DNA assemblies generated from long, uncorrected sequencing reads.

        !!! techdetails "Racon Technical Details"
            |  | Links |
            | --- | --- |
            | WDL Task | [task_racon.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/polishing/task_racon.wdl) |
            | Software Source Code | [Racon on GitHub](https://github.com/lbcb-sci/racon) |
            | Software Documentation | [Racon Documentation](https://github.com/lbcb-sci/racon#racon) |
            | Original Publication(s) | [Fast and accurate de novo genome assembly from long uncorrected reads](https://genome.cshlp.org/content/27/5/737) |

    ??? toggle "`Filter Contigs`: Filter contigs below a threshold length and remove homopolymer contigs"
        This task filters the created contigs based on a user-defined minimum length threshold (default of 1000) and eliminates homopolymer contigs (contigs of any length that consist of a single nucleotide). This ensures high-quality assemblies by retaining only contigs that meet specified criteria. Detailed metrics on contig counts and sequence lengths before and after filtering are provided in the output.

        !!! techdetails "Filter Contigs Technical Details" 
            |  | Links |
            | --- | --- |
            | WDL Task | [task_filter_contigs.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_filter_contigs.wdl) |
        
    ??? toggle "`Dnaapler`: Final Assembly Orientation"
        Dnaapler reorients contigs to start at specific reference points. Dnaapler supports the following modes, which can be indicated by filling the `dnaapler_mode` input variable with the desired mode. The default is `all`, which reorients contigs to start with `dnaA`, `terL`, `repA`, or `COG1474`.

        - **all**: Reorients contigs to start with `dnaA`, `terL`, `repA`, or `COG1474` (_Default_)
        - **chromosome**: Reorients to begin with the `dnaA` chromosomal replication initiator gene, commonly used for bacterial chromosome assemblies.
        - **plasmid**: Reorients to start with the `repA` plasmid replication initiation gene, ideal for plasmid assemblie
        - **phage**: Reorients to start with the `terL` large terminase subunit gene, used for bacteriophage assemblies
        - **archaea**: Reorients to start with the `COG1474` archaeal Orc1/cdc6 gene, relevant for archaeal assemblies
        - **custom**: Reorients based on a user-specified gene in amino acid FASTA format for experimental or unique workflows
        - **mystery**: Reorients to start with a random CDS for exploratory purposes
        - **largest**: Reorients to start with the largest CDS in the assembly, often useful for poorly annotated genomes
        - **nearest**: Reorients to start with the first CDS nearest to the sequence start, resolving CDS breakpoints
        - **bulk**: Processes multiple contigs to start with the desired start gene (`dnaA`, `terL`, `repA`, or custom)

        !!! techdetails "Dnaapler Technical Details"
            |  | Links |
            | --- | --- |
            | WDL Task | [task_dnaapler.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_dnaapler.wdl) |
            | Software Source Code | [Dnaapler on GitHub](https://github.com/gbouras13/dnaapler) |
            | Software Documentation | [Dnaapler Documentation](https://github.com/gbouras13/dnaapler?tab=readme-ov-file#dnaapler) |
            | Original Publication(s) | [Dnaapler: a tool to reorient circular microbial genomes](https://joss.theoj.org/papers/10.21105/joss.05968) |

#### Post-Assembly Tasks (performed for all taxa)

??? task  "`quast`: Assembly Quality Assessment"

    QUAST stands for QUality ASsessment Tool. It evaluates genome/metagenome assemblies by computing various metrics without a reference being necessary. It includes useful metrics such as number of contigs, length of the largest contig and N50. 

    !!! techdetails "QUAST Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_quast.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_quast.wdl) |
        | Software Source Code | [QUAST on GitHub](https://github.com/ablab/quast) |
        | Software Documentation | <https://quast.sourceforge.net/> |
        | Original Publication(s) | [QUAST: quality assessment tool for genome assemblies](https://academic.oup.com/bioinformatics/article/29/8/1072/228832) |

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

??? task "`MUMmer_ANI`: Average Nucleotide Identity (optional)"

    Average Nucleotide Identity (ANI) is a useful approach for taxonomic identification. The higher the percentage ANI of a query sequence to a given reference genome, the more likely the sequence is the same taxa as the reference. 

    ANI is calculated in TheiaProk using [a perl script written by Lee Katz](https://github.com/lskatz/ani-m) (ani-m.pl). This uses [MUMmer](http://mummer.sourceforge.net/) to rapidly align entire query assemblies to one or more reference genomes. By default, TheiaProk uses a set of 43 reference genomes in [RGDv2](https://github.com/StaPH-B/docker-builds/blob/master/build-files/fastani/1.34-RGDV2/RGDv2-metadata.tsv), a database containing genomes of enteric pathogens commonly sequenced by CDC EDLB & PulseNet participating laboratories. The user may also provide their own reference genome. After genome alignment with MUMmer, ani-m.pl calculates the average nucleotide identity and percent bases aligned between 2 genomes (query and reference genomes)

    The default database of reference genomes used is called "Reference Genome Database version 2" AKA "RGDv2". This database is composed of 43 enteric bacteria representing 32 species and is intended for identification of enteric pathogens and common contaminants. It contains six Campylobacter spp., three Escherichia/Shigella spp., one *Grimontia hollisae*, six *Listeria spp.*, one *Photobacterium damselae*, two *Salmonella spp.*, and thirteen *Vibrio spp.* 

    2 Thresholds are utilized to prevent false positive hits. The `ani_top_species_match` will only report a genus & species match if both thresholds are surpassed. Both of these thresholds are set to match those used in BioNumerics for PulseNet organisms.

    1. `ani_threshold` default value of 80.0
    2. `percent_bases_aligned_threshold` default value of 70.0

    For more information on RGDv2 database of reference genomes, please see [the publication here.](https://www.frontiersin.org/articles/10.3389/fmicb.2023.1225207/full)

    !!! techdetails "MUMmer_ANI Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_mummer_ani.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/advanced_metrics/task_mummer_ani.wdl) |
        | Software Source Code | [ani-m](https://github.com/lskatz/ani-m), [MUMmer](https://github.com/mummer4/mummer) |
        | Software Documentation | [ani-m](https://github.com/lskatz/ani-m), [MUMmer](https://mummer.sourceforge.net/) |
        | Original Publication(s) | [MUMmer4: A fast and versatile genome alignment system](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005944) |
        | Publication about RGDv2 database | https://www.frontiersin.org/articles/10.3389/fmicb.2023.1225207/full |

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

??? task "`KmerFinder`: Taxon Assignment (optional)"

    The `KmerFinder` method predicts prokaryotic species based on the number of overlapping (co-occurring) *k*-mers, i.e., 16-mers, between the query genome and genomes in a reference database.

    !!! techdetails "KmerFinder Technical Details"        
        
        |  | Links |
        | --- | --- |
        | Task | [task_kmerfinder.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/contamination/task_kmerfinder.wdl) |
        | Software Source Code | https://bitbucket.org/genomicepidemiology/kmerfinder |
        | Software Documentation | https://cge.food.dtu.dk/services/KmerFinder/instructions.php |
        | Original Publication(s) | [**Benchmarking of Methods for Genomic Taxonomy**](https://journals.asm.org/doi/full/10.1128/jcm.02981-13?rfr_dat=cr_pub++0pubmed&url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org) |

??? task "`AMRFinderPlus`: AMR Genotyping (default)"

    NCBI's [AMRFinderPlus](https://github.com/ncbi/amr/wiki) is the default antimicrobial resistance (AMR) detection tool used in TheiaProk. ResFinder may be used alternatively and if so, AMRFinderPlus is not run. 

    AMRFinderPlus identifies acquired antimicrobial resistance (AMR) genes, virulence genes, and stress genes.  Such AMR genes confer resistance to antibiotics, metals, biocides, heat, or acid. For some taxa (see [here](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus#--organism-option)), AMRFinderPlus will provide taxa-specific results including filtering out genes that are almost ubiquitous in the taxa (intrinsic genes) and identifying resistance-associated point mutations.  In TheiaProk, the taxon used by AMRFinderPlus is specified based on the `gambit_predicted_taxon` or a user-provided `expected_taxon`.

    You can check if a gene or point mutation is in the AMRFinderPlus database [here](https://www.ncbi.nlm.nih.gov/pathogens/refgene/#), find the sequences of reference genes [here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA313047), and search the query Hidden Markov Models (HMMs) used by AMRFinderPlus to identify AMR genes and some stress and virulence proteins ([here](https://www.ncbi.nlm.nih.gov/pathogens/hmm/)). The AMRFinderPlus database is updated frequently. You can ensure you are using the most up-to-date version by specifying the docker image as a workflow input. You might like to save this docker image as a workspace data element to make this easier.

    !!! techdetails "AMRFinderPlus Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_amrfinderplus.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/drug_resistance/task_amrfinderplus.wdl) |
        | Software Source Code | [amr on GitHub](https://github.com/ncbi/amr) |
        | Software Documentation | https://github.com/ncbi/amr/wiki |
        | Original Publication(s) | [AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8208984/) |

??? task "`ResFinder`: AMR Genotyping & Shigella XDR phenotype prediction (alternative)"

    The `ResFinder` task is an alternative to using AMRFinderPlus for detection and identification of AMR genes and resistance-associated mutations.

    This task runs the Centre for Genomic Epidemiology (CGE) ResFinder tool to identify acquired antimicrobial resistance. It can also run the CGE PointFinder tool if the `call_pointfinder` variable is set with to `true`. The databases underlying the task are different to those used by AMRFinderPlus.

    The default thresholds for calling AMR genes are 90% identity and 50% coverage of the reference genes (expressed as a fraction in workflow inputs: 0.9 & 0.5). These are the same thresholds utilized in BioNumerics for calling AMR genes.

    Organisms currently support by PointFinder for mutational-based predicted resistance:

    - Campylobacter coli & C. jejuni
    - Enterococcus faecalis
    - Enterococcus faecium
    - Escherichia coli & Shigella spp.
    - Helicobacter pylori
    - Neisseria gonorrhoeae
    - Klebsiella
    - Mycobacterium tuberculosis
    - Salmonella spp.
    - Staphylococcus aureus

    **XDR Shigella prediction**

    The `ResFinder` Task also has the ability to predict whether or not a sample meets the CDC's definition for extensively drug-resistant (XDR) Shigella. 

    > *CDC defines XDR Shigella bacteria as strains that are resistant to all commonly recommended empiric and alternative antibiotics — azithromycin, ciprofloxacin, ceftriaxone, trimethoprim-sulfamethoxazole (TMP-SMX), and ampicillin. [Link to CDC HAN](https://emergency.cdc.gov/han/2023/han00486.asp) where this definition is found.*
    
    A sample is required to meet **all 7 criteria** in order to be predicted as `XDR Shigella` 

    1. The GAMBIT task in the workflow must identify the sample as `Shigella` OR the user must input the word `Shigella` somewhere within the input String variable called `expected_taxon`. This requirement serves as the identification of a sample to be of the Shigella genus.
    2. Resfinder or PointFinder predicted resistance to **Ampicillin**
    3. Resfinder or PointFinder predicted resistance to **Azithromycin**
    4. Resfinder or PointFinder predicted resistance to **Ciprofloxacin**
    5. Resfinder or PointFinder predicted resistance to **Ceftriazone**
    6. Resfinder or PointFinder predicted resistance to **Trimethoprim**
    7. Resfinder or PointFinder predicted resistance to **Sulfamethoxazole**

    There are 3 potential outputs for the **`resfinder_predicted_xdr_shigella`** output string**:**

    - **`Not Shigella based on gambit_predicted_taxon or user input`**
    - **`Not XDR Shigella`** for samples identified as Shigella by GAMBIT or user input BUT does ResFinder did not predict resistance to **all 6 drugs in XDR definition**
    - **`XDR Shigella`** meaning the sample was identified as Shigella and ResFinder/PointFinder did predict resistance to ceftriazone, azithromycin, ciprofloxacin, trimethoprim, sulfamethoxazole, and ampicillin.
    
    !!! techdetails "ResFinder Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_resfinder.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/drug_resistance/task_resfinder.wdl) |
        | Software Source Code | https://bitbucket.org/genomicepidemiology/resfinder/src/master/ |
        | Software Documentation | https://bitbucket.org/genomicepidemiology/resfinder/src/master/ |
        | ResFinder database | https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/ |
        | PointFinder database | https://bitbucket.org/genomicepidemiology/pointfinder_db/src/master/ |
        | Web-server | https://cge.food.dtu.dk/services/ResFinder/ |
        | Original Publication(s) | [ResFinder 4.0 for predictions of phenotypes from genotypes](https://academic.oup.com/jac/article/75/12/3491/5890997) |

??? task "`TS_MLST`: MLST Profiling"

    [Multilocus sequence typing (MLST)](https://doi.org/10.1073/pnas.95.6.3140) is a typing method reflecting population structure. It was developed as a portable, unambiguous method for global epidemiology using PCR, but can be applied to whole-genome sequences *in silico*. MLST is commonly used for pathogen surveillance, ruling out transmission, and grouping related genomes for comparative analysis.

    MLST schemes are taxa-specific. Each scheme uses fragments of typically 7 housekeeping genes ("loci") and has a database associating an arbitrary number with each distinct allele of each locus. Each unique combination of alleles ("allelic profile") is assigned a numbered sequence type (ST). Significant diversification of genomes is captured by changes to the MLST loci via mutational events creating new alleles and STs, or recombinational events replacing the allele and changing the ST. Relationships between STs are based on the number of alleles they share. Clonal complexes share a scheme-specific number of alleles (usually for five of the seven loci).

    !!! tip "MLST Limitations"
        Some taxa have multiple MLST schemes, and some MLST schemes are insufficiently robust.

    TheiaProk uses [the MLST tool developed by Torsten Seeman](https://github.com/tseemann/mlst) to assess MLST using traditional [PubMLST](https://pubmlst.org/) typing schemes. 

    ??? toggle "Interpretation of MLST results"
        
        Each MLST results file returns the ST and allele results for one sample. If the alleles and ST are correctly assigned, only a single integer value will be present for each. If an ST cannot be assigned, multiple integers or additional characters will be shown, representing the issues with assignment as described [here](https://github.com/tseemann/mlst/tree/v2.22.0#missing-data).
        
    ??? toggle "Identifying novel alleles and STs"
        
        The MLST schemes used in TheiaProk are curated on the PubMLST website.If you identify novel alleles or allelic profiles in your data using TheiaProk's MLST task, you can get these assigned via PubMLST:
        
        1. Check that the novel allele or ST has not already been assigned a type on PubMLST. 
            1. Download the assembly file from Terra for your sample with the novel allele or ST
            2. Go to the [PubMLST webpage for the organism of interest](https://pubmlst.org/organisms) 
            3. Navigate to the organism "Typing" page 
            4. Under "Query a sequence" choose "Single sequence" (e.g., [this](https://pubmlst.org/bigsdb?db=pubmlst_hinfluenzae_seqdef&page=sequenceQuery) is the page for _H. influenzae_), select the MLST scheme under "Please select locus/scheme", upload the assembly fasta file, and click submit.
            5. Results will be returned lower on the page.
        2. If the allele or ST has not been typed previously on the PubMLST website (step 1), new allele or ST numbers can be assigned using instructions [here](https://pubmlst.org/submit-data).
        
    ??? toggle "Taxa with multiple MLST schemes"
        
        As default, the MLST tool automatically detects the genome's taxa to select the MLST scheme. 
        
        Some taxa have multiple MLST schemes, e.g. the *Escherichia* and Leptospira genera,  *Acinetobacter baumannii, Clostridium difficile* and *Streptococcus thermophilus.* Only one scheme will be used by default.
        
        Users may specify the scheme as an optional workflow input using the `scheme` variable of the `ts_mlst` task. Available schemes are listed [here](https://github.com/tseemann/mlst/blob/master/db/scheme_species_map.tab) and the scheme name should be provided in quotation marks ("...").
        
        If results from multiple MLST schemes are required for the same sample, TheiaProk can be run multiple times specifying non-default schemes. After the first run, output attributes for the workflow (i.e. output column names) must be amended to prevent results from being overwritten. Despite re-running the whole workflow, unmodified tasks will return cached outputs, preventing redundant computation.
        
    !!! techdetails "TS_MLST Technical Details"        
        |  | Links |
        | --- | --- |
        | Task | [task_ts_mlst.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/multi/task_ts_mlst.wdl) |
        | Software Source Code | [mlst](https://github.com/tseemann/mlst) |
        | Software Documentation | [mlst](https://github.com/tseemann/mlst) |

??? task "`Prokka`: Assembly Annotation (default)"

    Assembly annotation is available via `Prokka` as default, or alternatively via `Bakta`. When Prokka annotation is used, Bakta is not.

    [`Prokka`](https://github.com/tseemann/prokka) is a prokaryotic genome annotation tool used to identify and describe features of interest within the genome sequence. Prokka annotates there genome by querying databases described [here](https://github.com/tseemann/prokka#databases).

    !!! techdetails "Prokka Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_prokka.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/annotation/task_prokka.wdl) |
        | Software Source Code | [prokka](https://github.com/tseemann/prokka) |
        | Software Documentation | [prokka](https://github.com/tseemann/prokka) |
        | Original Publication(s) | [Prokka: rapid prokaryotic genome annotation](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517?login=false) |

??? task "`Bakta`: Assembly Annotation (alternative)"

    Assembly annotation is available via `Bakta` as an alternative to `Prokka`. When `Bakta` annotation is used, `Prokka` is not.

    `Bakta` is intended for annotation of Bacteria and plasmids only, and is best described [here](https://github.com/oschwengers/bakta#description)!

    In addition to the standard annotation outputs, `Bakta` also provides a plot summarizing the annotation results, which can be useful for visualizing genome features.

    **Bakta Database Options**

    `Bakta` supports three database configurations:

    **Light** Database: Optimized for faster performance and lower resource usage, with a focused set of core reference data for most bacterial genome annotations. Recommended for quick annotations or limited computational resources. Specify "light" for the `bakta_db` input.

    **Full** Database (default): Comprehensive with extensive reference annotations, suitable for detailed and accurate annotations. Specify "full" for the `bakta_db` input.

    **Custom** Database: Allows users to provide a Bakta-compatible database stored in Google Cloud Storage Must be a .tar.gz archive containing a properly formatted Bakta database with a valid version.json Follow the [Bakta database documentation](https://github.com/oschwengers/bakta#database) for detailed formatting requirements. Example: `"bakta_db": "gs://my-bucket/custom_bakta_db.tar.gz"`

    !!! techdetails "Bakta Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_bakta.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/annotation/task_bakta.wdl) |
        | Software Source Code | [bakta](https://github.com/oschwengers/bakta) |
        | Software Documentation | <https://github.com/oschwengers/bakta> |
        | Original Publication(s) | [Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000685) |

??? task "`PlasmidFinder`: Plasmid Identification"

    [`PlasmidFinder`](https://github.com/genomicepidemiology/plasmidfinder) detects plasmids in totally- or partially-sequenced genomes, and identifies the closest plasmid type in the database for typing purposes.

    ??? toggle "What are plasmids?"
        
        Plasmids are double-stranded circular or linear DNA molecules that are capable of replication independently of the chromosome and may be transferred between different species and clones. Many plasmids contain resistance or virulence genes, though some do not clearly confer an advantage to their host bacterium.
        
    !!! techdetails "PlasmidFinder Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_plasmidfinder.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/plasmid_detection/task_plasmidfinder.wdl) |
        | Software Source Code | https://bitbucket.org/genomicepidemiology/plasmidfinder/src/master/ |
        | Software Documentation | https://bitbucket.org/genomicepidemiology/plasmidfinder/src/master/ |
        | Original Publication(s) | [In Silico Detection and Typing of Plasmids using PlasmidFinder and Plasmid Multilocus Sequence Typing](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4068535/) |

??? task "**`qc_check`: Check QC Metrics Against User-Defined Thresholds (optional)**"

    The `qc_check` task compares generated QC metrics against user-defined thresholds for each metric. This task will run if the user provides a `qc_check_table` .tsv file. If all QC metrics meet the threshold, the `qc_check` output variable will read `QC_PASS`. Otherwise, the output will read `QC_NA` if the task could not proceed or `QC_ALERT` followed by a string indicating what metric failed.

    The `qc_check` task applies quality thresholds according to the sample taxa. The sample taxa is taken from the `gambit_predicted_taxon` value inferred by the GAMBIT module OR can be manually provided by the user using the `expected_taxon` workflow input.

    ??? toggle "Formatting the _qc_check_table.tsv_"

        - The first column of the qc_check_table lists the taxa that the task will assess and the header of this column must be "taxon".
        - Any genus or species can be included as a row of the qc_check_table. However, these taxa must ^^**uniquely**^^ match the sample taxa, meaning that the file can include multiple species from the same genus (Vibrio_cholerae and Vibrio_vulnificus), but not both a genus row and species within that genus (Vibrio and Vibrio cholerae). **The taxa should be formatted with the first letter capitalized and underscores in lieu of spaces.**
        - Each subsequent column indicates a QC metric and lists a threshold for each taxa that will be checked. **The column names must exactly match expected values, so we highly recommend copy and pasting from the template files below.**

    ??? toggle "Template _qc_check_table.tsv_ files"

        - TheiaProk_Illumina_PE: [theiaprok_illumina_pe_qc_check_template.tsv](../../assets/files/TheiaProk_Illumina_PE_qc_check_template.tsv)
        - TheiaProk_FASTA: [theiaprok_fasta_qc_check_template.tsv](../../assets/files/TheiaProk_FASTA_qc_check_template.tsv)

        !!! warning "Example Purposes Only"
            QC threshold values shown are for example purposes only and should not be presumed to be sufficient for every dataset.

    !!! techdetails "QC_Check Technical Details"    
        
        |  | Links |
        | --- | --- |
        | Task | [task_qc_check_phb.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/comparisons/task_qc_check_phb.wdl.wdl) |

??? task "`Taxon Tables`: Copy outputs to new data tables based on taxonomic assignment (optional)"

    The `taxon_tables` module, if enabled, will copy sample data to a different data table based on the taxonomic assignment. For example, if an *E. coli* sample is analyzed, the module will copy the sample data to a new table for *E. coli* samples or add the sample data to an existing table.

    !!! tip ""
        To implement the `taxon_tables` module, provide a file indicating data table names to copy samples of each taxa to in the `taxon_tables` input variable. No other input variables are needed.
        
        **Formatting the `taxon_tables` file**
        
        The `taxon_tables`  file must be uploaded a Google storage bucket that is accessible by Terra and should be in the format below. Briefly, the bacterial genera or species should be listed in the leftmost column with the name of the data table to copy samples of that taxon to in the rightmost column.
        
        | taxon | taxon_table |
        | --- | --- |
        | Listeria_monocytogenes | lmonocytogenes_specimen |
        | Salmonella | salmonella_specimen |
        | Escherichia | ecoli_specimen |
        | Shigella | shigella_specimen |
        | Streptococcus | strep_pneumo_specimen |
        | Legionella | legionella_specimen |
        | Klebsiella | klebsiella_specimen |
        | Mycobacterium | mycobacterium_specimen |
        | Acinetobacter | acinetobacter_specimen |
        | Pseudomonas | pseudomonas_specimen |
        | Staphylococcus | staphyloccus_specimen |
        | Neisseria | neisseria_specimen |
        
    !!! tip ""        
        There are no output columns for the taxon table task. The only output of the task is that additional data tables will appear for in the Terra workspace for samples matching a taxa in the `taxon_tables` file.

??? task "`Abricate`: Mass screening of contigs for antimicrobial and virulence genes (optional)"

    The `abricate` module, if enabled, will run abricate with the database defined in `abricate_db` to perform mass screening of contigs for antimicrobial resistance or virulence genes. It comes bundled with multiple databases: NCBI, CARD, ARG-ANNOT, Resfinder, MEGARES, EcOH, PlasmidFinder, Ecoli_VF and VFDB. It only detects acquired resistance genes, **NOT** point mutations

### Taxa-Specific Tasks

The TheiaProk workflows automatically activate taxa-specific sub-workflows after the identification of relevant taxa using `GAMBIT`. Alternatively, the user can provide the expected taxa in the `expected_taxon` workflow input to override the taxonomic assignment made by GAMBIT. Modules are launched for all TheiaProk workflows unless otherwise indicated.

??? toggle "_Acinetobacter baumannii_"
    ##### _Acinetobacter baumannii_ {#acinetobacter-baumannii}

    A number of approaches are available in TheiaProk for *A. baumannii* characterization.

    ??? task "`Kaptive`: Capsule and lipooligosaccharide outer core typing"
        
        The cell-surface capsular polysaccharide (CPS) of *Acinetobacter baumannii* can be used as an epidemiological marker. CPS varies in its composition and structure and is a key determinant in virulence and a target for non-antibiotic therapeutics. Specificity for non-antibiotic therapeutics (e.g. phage therapy) bear particular significance given the extent of antibiotic resistance found in this [ESKAPE](https://journals.asm.org/doi/10.1128/CMR.00181-19) pathogen. 
        
        Biosynthesis and export of CPS is encoded by genes clustering at the K locus (KL). Additional genes associated with CPS biosynthesis and export are sometimes found in other chromosomal locations. The full combination of these genes is summarized as a "[K type](https://www.biorxiv.org/content/10.1101/2022.05.19.492579v1)", described as a "predicted serotype associated with the best match locus". You can read more about this [here](https://github.com/katholt/Kaptive/wiki/Databases-distributed-with-Kaptive#acinetobacter-baunannii-k-and-oc-locus-databases).
        
        Previously, s[erotyping of *A. baumannii*](https://journals.asm.org/doi/10.1128/jcm.27.12.2713-2716.1989) focused on a major immunogenic polysaccharide which was considered the O antigen for the species. This serotyping approach appears to no longer be used and the serotyping [scheme has not been updated in over 20 years](https://www.karger.com/Article/Abstract/7300). Nonetheless, the O-antigen polysaccharide is attached to lipooligosaccharide, and the outer core (OC) of this lipooligosaccharide varies. Biosynthesis of the outer core lipooligosaccharide is encoded by a cluster of genes at the outer core (OC) locus.
        
        Variation in the KL and OCL can be characterized with the **Kaptive** tool and its associated [databases](https://github.com/katholt/Kaptive/wiki/Databases-distributed-with-Kaptive#acinetobacter-baunannii-k-and-oc-locus-databases) of numbered *A. baumannii* [K](https://github.com/katholt/Kaptive/blob/master/extras/Acinetobacter_baumannii_KL_reference_information.pdf) and [OC](https://github.com/katholt/Kaptive/blob/master/extras/Acinetobacter_baumannii_OCL_reference_information.pdf) locus variants. Kaptive takes in a genome assembly file (fasta), and assigns the K and OC locus to their numbered variants, provides K type and a description of genes in the K or OC loci and elsewhere in the chromosome, alongside metrics for quality of locus match. A description of [how Kaptive works](https://github.com/katholt/Kaptive/wiki/How-does-Kaptive-work%3F), [explanations of the full output reports](https://github.com/katholt/Kaptive/wiki/How-to-run#summary-table) which are provided in the Terra data table by TheiaProk and [resources for interpreting outputs](https://github.com/katholt/Kaptive/wiki/Interpreting-the-results) are available on the [Kaptive Wiki page](https://github.com/katholt/Kaptive/wiki/How-to-run#summary-table).
        
        !!! techdetails "Kaptive Technical Details"
            
            |  | Links |
            | --- | --- |
            | Task | [task_kaptive.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/acinetobacter/task_kaptive.wdl) |
            | Software Source Code | [Kaptive on GitHub](https://github.com/katholt/Kaptive/wiki) |
            | Software Documentation | https://github.com/katholt/Kaptive/wiki |
            | Orginal publications | [Identification of Acinetobacter baumannii loci for capsular polysaccharide (KL) and lipooligosaccharide outer core (OCL) synthesis in genome assemblies using curated reference databases compatible with Kaptive](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000339)<br>[An update to the database for Acinetobacter baumannii capsular polysaccharide locus typing extends the extensive and diverse repertoire of genes found at and outside the K locus](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000878) |
        
    ??? task "`AcinetobacterPlasmidTyping`: Acinetobacter plasmid detection"
        
        *Acinetobacter* plasmids are not included in the PlasmidFinder (see the relevant toggle in [this block](https://theiagen.github.io/public_health_bioinformatics/latest/workflows/genomic_characterization/theiaprok/#post-assembly-tasks-performed-for-all-taxa)) database. Instead, the [AcinetobacterPlasmidTyping](https://github.com/MehradHamidian/AcinetobacterPlasmidTyping) database contains variants of the plasmid *rep* gene for *A. baumannii* plasmid identification. When matched with >/= 95 % identity, this represents a typing scheme for *Acinetobacter baumannii* plasmids. In TheiaProk, we use the tool [abricate](https://github.com/tseemann/abricate) to query our assemblies against this database.
            
        The bioinformatics software for querying sample assemblies against the AcinetobacterPlasmidTyping database is [Abricate](https://github.com/tseemann/abricate). The WDL task simply runs abricate, and the Acinetobacter Plasmid database and default setting of 95% minimum identity are set in the [merlin magic sub-workflow](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_merlin_magic.wdl).

        !!! techdetails "AcinetobacterPlasmidTyping Technical Details"

            |  | Links |
            | --- | --- |
            | Task | [task_abricate.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/drug_resistance/task_abricate.wdl) |
            | Database and documentation | [https://github.com/MehradHamidian/AcinetobacterPlasmidTyping](https://github.com/MehradHamidian/AcinetobacterPlasmidTyping/tree/v1.0.0) |
            | Software Source Code and documentation | [abricate on GitHub](https://github.com/tseemann/abricate) |
            | Original Publication(s) | [Detection and Typing of Plasmids in *Acinetobacter baumannii* Using *rep* Genes Encoding Replication Initiation Proteins](https://journals.asm.org/doi/10.1128/spectrum.02478-22?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed) |
    
    ??? task "Acinetobacter MLST"
        
        Two MLST schemes are available for *Acinetobacter*. The Pasteur scheme is run by default, given [significant problems with the Oxford scheme have been described](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6510311/). Should users with to alternatively or additionally use the Oxford MLST scheme, see the section above on MLST. The Oxford scheme is activated in TheiaProk with the MLST `scheme` input as "abaumannii".
        
    ??? task "*bla*OXA-51-like gene detection"
        
        The *bla*OXA-51-like genes, also known as _oxaAB_, are considered intrinsic to _Acinetobacter baumannii_ but are not found in other *Acinetobacter* species. **Identification of a *bla*OXA-51-like gene is therefore considered to confirm the species' identity as _A. baumannii_.** 
        
        NCBI's AMRFinderPlus, which is implemented as a core module in TheiaProk, detects the *bla*OXA-51-like genes. This may be used to confirm the species, in addition to the GAMBIT taxon identification. The *bla*OXA-51-like genes act as carbapenemases when an IS*Aba1* is found 7 bp upstream of the gene. Detection of this IS is not currently undertaken in TheiaProk.

??? toggle "_Escherichia_ or _Shigella_ spp."
    ##### _Escherichia_ or _Shigella_ spp. {#escherichia-or-shigella}

    The *Escherichia* and *Shigella* genera are [difficult to differentiate as they do not comply with genomic definitions of genera and species](https://www.sciencedirect.com/science/article/abs/pii/S1286457902016374). Consequently, when either _Escherichia_ or _Shigella_ are identified by GAMBIT, all tools intended for these taxa are used. 

    `SerotypeFinder` and `ECTyper` are intended for analysis of *E. coli*. Both tools are used as there are occasional discrepancies between the serotypes predicted. This primarily arises due to differences in the databases used by each tool.

    ??? task "`SerotypeFinder`: Serotyping"
        
        [SerotypeFinder](https://bitbucket.org/genomicepidemiology/serotypefinder/src/master/), from the Centre for Genomic Epidemiology (CGE), identifies the serotype of total or partially-sequenced isolates of *E. coli*.
        
        !!! techdetails "SerotypeFinder Technical Details"
            
            |  | Links |
            | --- | --- |
            | Task | [task_serotypefinder.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/escherichia_shigella/task_serotypefinder.wdl) |
            | Software Source Code | https://bitbucket.org/genomicepidemiology/serotypefinder/src/master/ |
            | Software Documentation | https://bitbucket.org/genomicepidemiology/serotypefinder/src/master/ |
            | Original Publication(s) | [Rapid and Easy In Silico Serotyping of Escherichia coli Isolates by Use of Whole-Genome Sequencing Data](https://journals.asm.org/doi/10.1128/JCM.00008-15) |
        
    ??? task "`ECTyper`: Serotyping"
        
        [ECTyper](https://github.com/phac-nml/ecoli_serotyping) is a serotyping module for *E. coli*. In TheiaProk, we are using assembly files as input.
        
        !!! techdetails "ECTyper Technical Details"
            
            |  | Links |
            | --- | --- |
            | Task | [task_ectyper.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/escherichia_shigella/task_ectyper.wdl) |
            | Software Source Code | [ECTyper on GitHub](https://github.com/phac-nml/ecoli_serotyping) |
            | Software Documentation | [ECTyper on GitHub](https://github.com/phac-nml/ecoli_serotyping) |
            | Orginal publication | [ECTyper: in silico Escherichia coli serotype and species prediction from raw and assembled whole-genome sequence data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8767331/) |

    `VirulenceFinder` identifies virulence genes in total or partial sequenced isolates of bacteria. Currently, only *E. coli* is supported in TheiaProk workflows. 

    ??? task "`VirulenceFinder`: Virulence gene identification"
        
        VirulenceFinder in TheiaProk is only run on assembly files due to issues regarding discordant results when using read files on the web application versus the command-line.
        
        !!! techdetails "VirulenceFinder Technical Details"
            
            |  | Links |
            | --- | --- |
            | Task | [task_virulencefinder.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/escherichia_shigella/task_virulencefinder.wdl) |
            | Software Source Code | [**VirulenceFinder**](https://bitbucket.org/genomicepidemiology/virulencefinder/src/master/) |
            | Software Documentation | [**VirulenceFinder**](https://bitbucket.org/genomicepidemiology/virulencefinder/src/master/) |
            | Original Publication(s) | [Real-time whole-genome sequencing for routine typing, surveillance, and outbreak detection of verotoxigenic Escherichia co](https://pubmed.ncbi.nlm.nih.gov/24574290/) |

    `ShigaTyper` and `ShigEiFinder` are intended for differentiation and serotype prediction for any *Shigella* species and Enteroinvasive *Escherichia coli* (EIEC). You can read about differences between these [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC517479/) and [here](https://www.microbiologyresearch.org/content/journal/micro/10.1099/00221287-144-9-2667). ShigEiFinder can be run using either the assembly (default) or reads. These tasks will report if the samples are neither *Shigella* nor EIEC.

    ??? task "`ShigaTyper`: *Shigella*/EIEC differentiation and serotyping ==_for Illumina and ONT only_=="
        
        ShigaTyper predicts *Shigella* spp. serotypes from Illumina or ONT read data. If the genome is not *Shigella* or EIEC, the results from this tool will state this. In the notes it provides, it also reports on the presence of *ipaB* which is suggestive of the presence of the "virulent invasion plasmid".
        
        !!! techdetails "ShigaTyper Technical Details"
            
            |  | Links |
            | --- | --- |
            | Task | [task_shigatyper.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/escherichia_shigella/task_shigatyper.wdl) |
            | Software Source Code | [ShigaTyper on GitHub](https://github.com/CFSAN-Biostatistics/shigatyper) |
            | Software Documentation | https://github.com/CFSAN-Biostatistics/shigatyper |
            | Origin publication | [In Silico Serotyping Based on Whole-Genome Sequencing Improves the Accuracy of Shigella Identification](https://doi.org/10.1128/AEM.00165-19) |
        
    ??? task "`ShigEiFinder`: *Shigella*/EIEC differentiation and serotyping ==_using the assembly file as input_=="
        
        ShigEiFinder differentiates *Shigella* and enteroinvasive *E. coli* (EIEC) using cluster-specific genes, identifies some serotypes based on the presence of O-antigen and H-antigen genes, and predicts the number of virulence plasmids. The `shigeifinder` task operates on assembly files.
        
        !!! techdetails "ShigEiFinder Technical Details"
            
            |  | Links |
            | --- | --- |
            | Task | [task_shigeifinder.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/escherichia_shigella/task_shigeifinder.wdl) |
            | Software Source Code | [ShigEiFinder on GitHub](https://github.com/LanLab/ShigEiFinder) |
            | Software Documentation | [ShigEiFinder on GitHub](https://github.com/LanLab/ShigEiFinder) |
            | Origin publication | [Cluster-specific gene markers enhance Shigella and enteroinvasive Escherichia coli in silico serotyping](https://pubmed.ncbi.nlm.nih.gov/34889728/) |
    
    ??? task "`ShigEiFinder_reads`: *Shigella*/EIEC differentiation and serotyping using Illumina read files as input (optional) ==_ for Illumina data only_=="

        ShigEiFinder differentiates *Shigella* and enteroinvasive *E. coli* (EIEC) using cluster-specific genes, identifies some serotypes based on the presence of O-antigen and H-antigen genes, and predicts the number of virulence plasmids. The `shigeifinder_reads` task performs on read files.
        
        !!! techdetails "ShigEiFinder_reads Technical Details"
            |  | Links |
            | --- | --- |
            | Task | [task_shigeifinder.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/escherichia_shigella/task_shigeifinder.wdl#L81) |
            | Software Source Code | [ShigEiFinder on GitHub](https://github.com/LanLab/ShigEiFinder) |
            | Software Documentation | [ShigEiFinder on GitHub](https://github.com/LanLab/ShigEiFinder) |
            | Origin publication | [Cluster-specific gene markers enhance Shigella and enteroinvasive Escherichia coli in silico serotyping](https://pubmed.ncbi.nlm.nih.gov/34889728/) |

    `SonneiTyper` is run only when GAMBIT predicts the *S. sonnei* species. This is the most common *Shigella* species in the United States.

    ??? task "`SonneiTyper`**: *Shigella sonnei* identification, genotyping, and resistance mutation identification ==_for Illumina and ONT data only_=="
        
        SonneiTyper identifies *Shigella sonnei,* and uses **single-nucleotide variants for genotyping and prediction of quinolone resistance in *gyrA* (S83L, D87G, D87Y) and *parC* (S80I). Outputs are provided in [this](https://github.com/katholt/sonneityping#example-output) format.
        
        SonneiTyper is a wrapper script around another tool, Mykrobe, that analyses the *S. sonnei* genomes.

        !!! techdetails "SonneiTyper Technical Details"

            |  | Links |
            | --- | --- |
            | Task | [task_sonneityping.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/escherichia_shigella/task_sonneityping.wdl) |
            | Software Source Code | [Mykrobe](https://github.com/Mykrobe-tools/mykrobe), [sonneityping](https://github.com/katholt/sonneityping) |
            | Software Documentation | https://github.com/Mykrobe-tools/mykrobe/wiki, [sonneityping](https://github.com/katholt/sonneityping) |
            | Original Publication(s) | [Global population structure and genotyping framework for genomic surveillance of the major dysentery pathogen, *Shigella sonnei*](https://www.nature.com/articles/s41467-021-22700-4) |

    **Shigella XDR prediction.** Please see the documentation section above for ResFinder for details regarding this taxa-specific analysis. 

    ??? task "`StxTyper`: Identification and typing of Shiga toxin (Stx) genes ==_using the assembly file as input_=="
        
        StxTyper screens bacterial genome assemblies for shiga toxin genes and subtypes them into known subtypes and also looks for novel subtypes in cases where the detected sequences diverge from the reference sequences.
        
        Shiga toxin is the main virulence factor of Shiga-toxin-producing E. coli (STEC), though these genes are also found in Shigella species as well as some other genera more rarely, such as Klebsiella. [Please see this review paper that describes shiga toxins in great detail.](https://doi.org/10.3390/microorganisms12040687)

        !!! tip "Running StxTyper via the TheiaProk workflows"
            The TheiaProk workflow will automatically run `stxtyper` on all E. coli and Shigella spp. samples, but ==*the user can opt-in to running the tool on any sample by setting the optional input variable `call_stxtyper` to `true` when configuring the workflow.*==
        
        Generally, `stxtyper` looks for _stxA_ and _stxB_ subunits that compose a complete operon. The A subunit is longer (in amino acid length) than the B subunit. Stxtyper attempts to detect these, compare them to a database of known sequences, and type them based on amino acid composition.  There typing algorithm and rules defining how to type these genes & operons will be described more completely in a publication that will be available in the future.
        
        The `stxtyper_report` output TSV is provided in [this output format.](https://github.com/ncbi/stxtyper/tree/v1.0.24?tab=readme-ov-file#output)

        This tool has been incorporated into v4.0.3 of AMRFinderPlus and runs behind-the-scenes when the user (or in this case, the TheiaProk workflow) provides the `amrfinder --organism Escherichia --plus` options.

        !!! techdetails "StxTyper Technical Details"

            |  | Links |
            | --- | --- |
            | Task | [task_stxtyper.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/escherichia_shigella/task_stxtyper.wdl) |
            | Software Source Code | [ncbi/stxtyper GitHub repository](https://github.com/ncbi/stxtyper) |
            | Software Documentation | [ncbi/stxtyper GitHub repository](https://github.com/ncbi/stxtyper) |
            | Original Publication(s) | No publication currently available, as this is a new tool. One will be available in the future. |

??? toggle "_Haemophilus influenzae_"
    ##### _Haemophilus influenzae_ {#haemophilus-influenzae}
    ??? task "`hicap`: Sequence typing"
        Identification of _cap_ locus serotype in _Haemophilus influenzae_ assemblies with [hicap](https://github.com/scwatts/hicap).

        The _cap_ locus of _H. influenzae_ is categorised into 6 different groups based on serology (a-f). There are three functionally distinct regions of the _cap_ locus, designated `region I`, `region II`, and `region III`. Genes within `region I` (`bexABCD`) and `region III` (`hcsAB`) are associated with transport and post-translation modification. The `region II` genes encode serotype-specific proteins, with each serotype (a-f) having a distinct set of genes. _cap_ loci are often subject to structural changes (e.g. duplication, deletion) making the process of *in silico* typing and characterisation of loci difficult.
        
        `hicap` automates the identification of the _cap_ locus, describes the structural layout, and performs *in silico* serotyping.
        
        !!! techdetails "hicap Technical Details"
            
            |  | Links |
            | --- | --- |
            | Task | [task_hicap.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/haemophilus/task_hicap.wdl) |
            | Software Source Code | [hicap on GitHub](https://github.com/scwatts/hicap) |
            | Software Documentation | [hicap on GitHub](https://github.com/scwatts/hicap) |
            | Original Publication(s) | [hicap: In Silico Serotyping of the Haemophilus influenzae Capsule Locus](https://doi.org/10.7717/peerj.5261) |

??? toggle "_Klebsiella_ spp."
    ##### _Klebsiella_ spp. {#klebsiella}
    ??? task "`Kleborate`: Species identification, MLST, serotyping, AMR and virulence characterization"

        [Kleborate](https://github.com/katholt/Kleborate) is a tool to identify the *Klebsiella* species, MLST sequence type, serotype, virulence factors (ICE_Kp_ and plasmid associated), and AMR genes and mutations. Serotyping is based on the capsular (K antigen) and lipopolysaccharide (LPS) (O antigen) genes. The resistance genes identified by Kleborate are described [here](https://github.com/katholt/Kleborate/wiki/Antimicrobial-resistance).
        
        !!! techdetails "Kleborate Technical Details"
            
            |  | Links |
            | --- | --- |
            | Task | [task_kleborate.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/klebsiella/task_kleborate.wdl) |
            | Software Source Code | [kleborate on GitHub](https://github.com/katholt/Kleborate) |
            | Software Documentation | https://github.com/katholt/Kleborate/wiki |
            | Orginal publication | [A genomic surveillance framework and genotyping tool for Klebsiella pneumoniae and its related species complex](https://www.nature.com/articles/s41467-021-24448-3)<br>[Identification of Klebsiella capsule synthesis loci from whole genome data](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000102) |

??? toggle "_Legionella pneumophila_"
    ##### _Legionella pneumophila_ {#legionella-pneumophila}
    ??? task "`Legsta`: Sequence-based typing"

        [Legsta](https://github.com/tseemann/legsta) performs a sequence-based typing of *Legionella pneumophila*, with the intention of being used for outbreak investigations.
        
        !!! techdetails "Legsta Technical Details"            
            
            |  | Links |
            | --- | --- |
            | Task | [task_legsta.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/legionella/task_legsta.wdl) |
            | Software Source Code | [Legsta](https://github.com/tseemann/legsta) |
            | Software Documentation | [Legsta](https://github.com/tseemann/legsta) |

??? toggle "_Listeria monocytogenes_"
    ##### _Listeria monocytogenes_ {#listeria-monocytogenes}
    ??? task "`LisSero`: Serogroup prediction"

        [LisSero](https://github.com/MDU-PHL/LisSero) performs serogroup prediction (1/2a, 1/2b, 1/2c, or 4b) for _Listeria monocytogenes_ based on the presence or absence of five genes, _lmo1118_, _lmo0737_, ORF2110, ORF2819, and _prs_. These do not predict somatic (O) or flagellar (H) biosynthesis.
        
        !!! techdetails "LisSero Technical Details"
            
            |  | Links |
            | --- | --- |
            | Task | [task_lissero.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/listeria/task_lissero.wdl) |
            | Software Source Code | [LisSero](https://github.com/MDU-PHL/LisSero) |
            | Software Documentation | [LisSero](https://github.com/MDU-PHL/LisSero) |

??? toggle "_Mycobacterium tuberculosis_"
    ##### _Mycobacterium tuberculosis_ {#mycobacterium-tuberculosis}
    ??? task "`TBProfiler`: Lineage and drug susceptibility prediction ==_for Illumina and ONT only_=="

        [TBProfiler](https://github.com/jodyphelan/TBProfiler) identifies *Mycobacterium tuberculosis* complex species, lineages, sub-lineages and drug resistance-associated mutations.
        
        !!! techdetails "TBProfiler Technical Details"
            
            |  | Links |
            | --- | --- |
            | Task | [task_tbprofiler.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/mycobacterium/task_tbprofiler.wdl) |
            | Software Source Code | [TBProfiler on GitHub](https://github.com/jodyphelan/TBProfiler) |
            | Software Documentation | https://jodyphelan.gitbook.io/tb-profiler/ |
            | Original Publication(s) | [Integrating informatics tools and portable sequencing technology for rapid detection of resistance to anti-tuberculous drugs](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0650-x) |
    
    ??? task "`tbp-parser`: Interpretation and Parsing of TBProfiler JSON outputs; ==_requires TBProfiler and `tbprofiler_additonal_outputs = true`_=="
    
        [tbp-parser](https://github.com/theiagen/tbp-parser/) adds useful drug resistance interpretation by applying expert rules and organizing the outputs from TBProfiler. Please note that this tool has **not** been tested on ONT data and although it is available, result accuracy should be considered carefully. To understand this module and its functions, [please examine the README found with the source code here](https://github.com/theiagen/tbp-parser/).
        
        !!! techdetails "tbp-parser Technical Details"
            
            |  | Links |
            | --- | --- |
            | Task | [task_tbp_parser.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/mycobacterium/task_tbp_parser.wdl) |
            | Software Source Code | [tbp-parser](https://github.com/theiagen/tbp-parser/) |
            | Software Documentation | [tbp-parser](https://theiagen.github.io/tbp-parser)  |

    ??? task "`Clockwork`: Decontamination of input read files ==_for Illumina PE only_=="
    
        [Clockwork](https://github.com/iqbal-lab-org/clockwork/wiki) decontaminates paired-end data by removing all reads that do not match the H37Rv genome or are unmapped.
        
        !!! techdetails "Clockwork Technical Details"

            |  | Links |
            | --- | --- |
            | Task | [task_clockwork.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/mycobacterium/task_clockwork.wdl) |
            | Software Source Code | [clockwork](https://github.com/iqbal-lab-org/clockwork) |
            | Software Documentation | <https://github.com/iqbal-lab-org/clockwork/wiki> |

??? toggle "_Neisseria_ spp."
    ##### _Neisseria_ spp. {#neisseria}

    ??? task "`amr_search`: _Neisseria gonorrhoeae_ antimicrobial resistance profiling"

        This task performs *in silico* antimicrobial resistance (AMR) profiling for *Neisseria gonorrhoeae* using **AMRsearch**, the primary tool used by [Pathogenwatch](https://pathogen.watch/) to genotype and infer antimicrobial resistance (AMR) phenotypes from assembled microbial genomes.

        **AMRsearch** screens against Pathogenwatch's library of curated genotypes and inferred phenotypes, developed in collaboration with community experts. Resistance phenotypes are determined based on both **resistance genes** and **mutations**, and the system accounts for interactions between multiple SNPs, genes, and suppressors. Predictions follow **S/I/R classification** (*Sensitive, Intermediate, Resistant*).

        The AMR search is conducted when *Neisseria gonorrhoeae* is identified as the taxon in *TheiaProk* workflows. The default database for *N. gonorrhoeae* is **485**.

        **Outputs:**

        - **JSON Output:** Contains the complete AMR profile, including detailed **resistance state**, detected **resistance genes/mutations**, and supporting **BLAST results**.

        - **CSV & PNG Tables:* Results are formatted into a **CSV file** and **PNG summary table** for easier visualization.

        !!! techdetails "amr_search Technical Details"    

            |  | Links |
            | --- | --- |
            | Task | [task_amr_search.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/drug_resistance/task_amr_search.wdl) |
            | Software Source Code | [AMRsearch](https://github.com/pathogenwatch-oss/amr-search) |
            | Software Documentation | [Pathogenwatch](https://cgps.gitbook.io/pathogenwatch) |
            | Original Publication(s) | [PAARSNP: *rapid genotypic resistance prediction for *Neisseria gonorrhoeae*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7545138/) |

    ??? task "`ngmaster`: _Neisseria gonorrhoeae_ sequence typing"

        NG-MAST is currently the most widely used method for epidemiological surveillance of *Neisseria gonorrhoea.* This tool is targeted at clinical and research microbiology laboratories that have performed WGS of *N. gonorrhoeae* isolates and wish to understand the molecular context of their data in comparison to previously published epidemiological studies. As WGS becomes more routinely performed, *NGMASTER*
         has been developed to completely replace PCR-based NG-MAST, reducing time and labour costs. 
        
        The NG-STAR offers a standardized method of classifying seven well-characterized genes associated antimicrobial resistance in *N. gonorrhoeae* (*penA, mtrR, porB, ponA, gyrA, parC* and 23S rRNA) to three classes of antibiotics (cephalosporins, macrolides and fluoroquinolones).
        
        ngmaster combines two tools: NG-MAST (*in silico* multi-antigen sequencing typing) and NG-STAR (sequencing typing for antimicrobial resistance).
        
        !!! techdetails "ngmaster Technical Details"    
            
            |  | Links |
            | --- | --- |
            | Task | [task_ngmaster.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/neisseria/task_ngmaster.wdl) |
            | Software Source Code | [ngmaster](https://github.com/MDU-PHL/ngmaster) |
            | Software Documentation | [ngmaster](https://github.com/MDU-PHL/ngmaster) |
            | Original Publication(s) | [NGMASTER: *in silico* multi-antigen sequence typing for *Neisseria gonorrhoeae*](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000076) |
        
    ??? task "`meningotype`: _Neisseria meningitidis_ serotyping"
        
        This tool performs serotyping, MLST, finetyping (of *porA*, *fetA*, and *porB*), and Bexsero Antigen Sequencing Typing (BAST). 
        
        !!! techdetails "meningotype Technical Details"
            
            |  | Links |
            | --- | --- |
            | Task | [task_meningotype.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/neisseria/task_meningotype.wdl) |
            | Software Source Code | [meningotype](https://github.com/MDU-PHL/meningotype) |
            | Software Documentation | [meningotype](https://github.com/MDU-PHL/meningotype) |

??? toggle "_Pseudomonas aeruginosa_"
    ##### _Pseudomonas aeruginosa_ {#pseudomonas-aeruginosa}

    ??? task "`pasty`: Serotyping"
        
        `pasty` is a tool for *in silico* serogrouping of *Pseudomonas aeruginosa* isolates. pasty was developed by Robert Petit, based on the [PAst](https://github.com/Sandramses/PAst) tool from the Centre for Genomic Epidemiology.
        
        !!! techdetails "pasty Technical Details"
            
            |  | Links |
            | --- | --- |
            | Task | [task_pasty.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/pseudomonas/task_pasty.wdl) |
            | Software Source Code | [pasty](https://github.com/rpetit3/pasty) |
            | Software Documentation | [pasty](https://github.com/rpetit3/pasty) |
            | Original Publication(s) | [Application of Whole-Genome Sequencing Data for O-Specific Antigen Analysis and In Silico Serotyping of Pseudomonas aeruginosa Isolates.](https://journals.asm.org/doi/10.1128/JCM.00349-16) |

??? toggle "_Salmonella_ spp."
    ##### _Salmonella_ spp. {#salmonella}

    Both SISTR and SeqSero2 are used for serotyping all *Salmonella* spp. Occasionally, the predicted serotypes may differ between SISTR and SeqSero2. When this occurs, differences are typically small and analogous, and are likely as a result of differing source databases. More information about Salmonella serovar nomenclature can be found [here](https://www.happykhan.com/posts/binfie-guide-serovar/). For *Salmonella* Typhi, genotyphi is additionally run for further typing.

    ??? task "`SISTR`: Salmonella serovar prediction"
        
        [SISTR](https://github.com/phac-nml/sistr_cmd) performs *Salmonella spp.* serotype prediction using antigen gene and cgMLST gene alleles. In TheiaProk. SISTR is run on genome assemblies, and uses the default database setting (smaller "centroid" alleles or representative alleles instead of the full set of cgMLST alleles). It also runs a QC mode to determine the level of confidence in the serovar prediction (see [here](https://github.com/phac-nml/sistr_cmd#qc-by-sistr_cmd---qc)).
        
        !!! techdetails "SISTR Technical Details"
            
            |  | Links |
            | --- | --- |
            | Task | [task_sistr.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/salmonella/task_sistr.wdl) |
            | Software Source Code | [SISTR](https://github.com/phac-nml/sistr_cmd) |
            | Software Documentation | [SISTR](https://github.com/phac-nml/sistr_cmd) |
            | Original Publication(s) | [The Salmonella In Silico Typing Resource (SISTR): an open web-accessible tool for rapidly typing and subtyping draft Salmonella genome assemblies.](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0147101) |

    ??? task "`SeqSero2`: Serotyping"
        
        [SeqSero2](https://github.com/denglab/SeqSero2) is a tool for *Salmonella* serotype prediction. In the TheiaProk Illumina and ONT workflows, SeqSero2 takes in raw sequencing reads and performs targeted assembly of serotype determinant alleles, which can be used to predict serotypes including contamination between serotypes. Optionally, SeqSero2 can take the genome assembly as input.
        
        !!! techdetails "SeqSero2 Technical Details"
            |  | Links |
            | --- | --- |
            | Task | [task_seqsero2.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/salmonella/task_seqsero2.wdl) |
            | Software Source Code | [SeqSero2](https://github.com/denglab/SeqSero2) |
            | Software Documentation | [SeqSero2](https://github.com/denglab/SeqSero2) |
            | Original Publication(s) | [Salmonella serotype determination utilizing high-throughput genome sequencing data.](https://journals.asm.org/doi/10.1128/JCM.00323-15)<br>[SeqSero2: rapid and improved Salmonella serotype determination using whole genome sequencing data.](https://journals.asm.org/doi/10.1128/AEM.01746-19) |
    
    ??? task "`genotyphi`: _Salmonella_ Typhi lineage, clade, subclade and plasmid typing, AMR prediction ==_for Illumina and ONT only_=="
        
        [`genotyphi`](https://github.com/katholt/genotyphi) is activated upon identification of the "Typhi" serotype by SISTR or SeqSero2. `genotyphi` divides the *Salmonella enterica* serovar Typhi population into detailed lineages, clades, and subclades. It also detects mutations in the quinolone-resistance determining regions, acquired antimicrobial resistance genes, plasmid replicons, and subtypes of the IncHI1 plasmid which is associated with multidrug resistance. 
        
        TheiaProk uses the [Mykrobe implementation](https://github.com/katholt/genotyphi/blob/main/README.md#mykrobe-implementation) of genotyphi that takes raw sequencing reads as input. 
        
        !!! techdetails "genotyphi Technical Details"
            
            |  | Links |
            | --- | --- |
            | Task | [task_genotyphi.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/salmonella/task_genotyphi.wdl) |
            | Software Source Code | [genotyphi](https://github.com/katholt/genotyphi) |
            | Software Documentation | https://github.com/katholt/genotyphi/blob/main/README.md#mykrobe-implementation |
            | Orginal publication(s) | [An extended genotyping framework for Salmonella enterica serovar Typhi, the cause of human typhoid](https://www.nature.com/articles/ncomms12827/)<br>[Five Years of GenoTyphi: Updates to the Global Salmonella Typhi Genotyping Framework](https://academic.oup.com/jid/article/224/Supplement_7/S775/6358992?login=false) |

??? toggle "_Staphyloccocus aureus_"
    ##### _Staphyloccocus aureus_ {#staphyloccocus-aureus}

    ??? task "`spatyper`: Sequence typing"
        
        Given a fasta file or multiple fasta files, this script identifies the repeats and the order and generates a *spa* type. The repeat sequences and repeat orders found on http://spaserver2.ridom.de/ are used to identify the spa type of each enriched sequence. Ridom *spa* type and the genomics repeat sequence are then reported back to the user.
        
        !!! techdetails "spatyper Technical Details"

            |  | Links |
            | --- | --- |
            | Task | [task_spatyper.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/staphylococcus/task_spatyper.wdl) |
            | Software Source Code | [spatyper](https://github.com/HCGB-IGTP/spaTyper) |
            | Software Documentation | [spatyper](https://github.com/HCGB-IGTP/spaTyper) |
            
    ??? task "`staphopia-sccmec`: Sequence typing"
        
        This tool assigns a SCCmec type by BLAST the SCCmec primers against an assembly. `staphopia-sccmec`reports `True` for exact primer matches and `False` for at least 1 base pair difference. The [Hamming Distance](https://en.wikipedia.org/wiki/Hamming_distance) is also reported.
        
        !!! techdetails "staphopia-sccmec Technical Details"           
            
            |  | Links |
            | --- | --- |
            | Task | [task_staphopiasccmec.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/staphylococcus/task_staphopiasccmec.wdl) |
            | Software Source Code | [staphopia-sccmec](https://github.com/staphopia/staphopia-sccmec) |
            | Software Documentation | [staphopia-sccmec](https://github.com/staphopia/staphopia-sccmec) |
            | Original Publication(s) | [*Staphylococcus aureus* viewed from the perspective of 40,000+ genomes](https://doi.org/10.7717/peerj.5261) |

    ??? task "`agrvate`: Sequence typing"
        
        This tool identifies the *agr* locus type and reports possible variants in the *agr* operon. AgrVATE accepts a *S. aureus* genome assembly as input and performs a kmer search using an Agr-group specific kmer database to assign the Agr-group. The *agr* operon is then extracted using *in-silico* PCR and variants are called using an Agr-group specific reference operon.
        
        !!! techdetails "agrvate Technical Details"
            
            |  | Links |
            | --- | --- |
            | Task | [task_agrvate.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/staphylococcus/task_agrvate.wdl) |
            | Software Source Code | [agrVATE](https://github.com/VishnuRaghuram94/AgrVATE) |
            | Software Documentation | [agrVATE](https://github.com/VishnuRaghuram94/AgrVATE) |
            | Original Publication(s) | [Species-Wide Phylogenomics of the *Staphylococcus aureus Agr* Operon Revealed Convergent Evolution of Frameshift Mutations](https://doi.org/10.1128/spectrum.01334-21) |

??? toggle "_Streptococcus pneumoniae_"
    ##### _Streptococcus pneumoniae_ {#streptococcus-pneumoniae}

    ??? task "`PopPUNK`: Global Pneumococcal Sequence Cluster typing"
        
        Global Pneumococcal Sequence Clusters (GPSC) define and name pneumococcal strains. GPSC designation is undertaken using the PopPUNK software and GPSC database as described in the file below, obtained from [here](https://www.pneumogen.net/gps/#/training#command-line).

        :file: [GPSC_README_PopPUNK2.txt](../../assets/files/GPSC_README_PopPUNK2.txt)
        
        !!! tip "Interpreting GPSC results"
            - In the `*_external_clusters.csv` novel clusters are assigned NA. For isolates that are assigned a novel cluster and pass QC, you can email [globalpneumoseq@gmail.com](mailto:globalpneumoseq@gmail.com) to have these novel clusters added to the database.
            - Unsampled diversity in the pneumococcal population may represent missing variation that links two GPS clusters. When this is discovered, GPSCs are merged and the merge history is indicated. For example, if GPSC23 and GPSC362 merged, the GPSC would be reported as GPSC23, with a merge history of GPSC23;362.
        
        !!! techdetails "PopPUNK Technical Details"
            
            |  | Links |
            | --- | --- |
            | Task | [task_poppunk_streppneumo.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/streptococcus/task_poppunk_streppneumo.wdl) |
            | GPSC database | <https://www.pneumogen.net/gps/#/training#command-line> |
            | Software Source Code | [PopPunk](https://github.com/bacpop/PopPUNK) |
            | Software Documentation | <https://poppunk.readthedocs.io/en/latest/> |
            | Original Publication(s) | [Fast and flexible bacterial genomic epidemiology with PopPUNK](https://genome.cshlp.org/content/29/2/304) |
        
    ??? task "`SeroBA`: Serotyping ==_for Illumina_PE only_=="
        
        Streptococcus pneumoniae serotyping is performed with SeroBA.
        
        !!! techdetails "SeroBA Technical Details"
            
            |  | Links |
            | --- | --- |
            | Task | [task_seroba.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/streptococcus/task_seroba.wdl) |
            | Software Source Code | [SeroBA](https://github.com/sanger-pathogens/seroba) |
            | Software Documentation | https://sanger-pathogens.github.io/seroba/ |
            | Original Publication(s) | [SeroBA: rapid high-throughput serotyping of Streptococcus pneumoniae from whole genome sequence data](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000186) |
    
    ??? task "`pbptyper`: Penicillin-binding protein genotyping"
       
        The Penicillin-binding proteins (PBP) are responsible for the minimum inhibitory concentration phenotype for beta-lactam antibiotic. In *Streptococcus pneumoniae*, these PBP genes can be identified and typed with PBPTyper. 
        
        !!! techdetails "pbptyper Technical Details"
            
            |  | Links |
            | --- | --- |
            | Task | [task_pbptyper.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/streptococcus/task_pbptyper.wdl) |
            | Software Source Code | [pbptyper](https://github.com/rpetit3/pbptyper) |
            | Software Documentation | [pbptyper](https://github.com/rpetit3/pbptyper) |
            | Original Publication(s) | [Penicillin-binding protein transpeptidase signatures for tracking and predicting β-lactam resistance levels in Streptococcus pneumoniae](https://journals.asm.org/doi/full/10.1128/mBio.00756-16) |

??? toggle "_Streptococcus pyogenes_"
    ##### _Streptococcus pyogenes_ {#streptococcus-pyogenes}
    ??? task "`emm-typing-tool`: Sequence typing ==_for Illumina_PE only_=="

        emm-typing of *Streptococcus pyogenes* raw reads. Assign emm type and subtype by querying the CDC M-type specific database. 
        
        !!! techdetails "emm-typing-tool Technical Details"            
            |  | Links |
            | --- | --- |
            | Task | [task_emmtypingtool.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/streptococcus/task_emmtypingtool.wdl) |
            | Software Source Code | [emm-typing-tool](https://github.com/ukhsa-collaboration/emm-typing-tool) |
            | Software Documentation | [emm-typing-tool](https://github.com/ukhsa-collaboration/emm-typing-tool) |

??? toggle "_Vibrio_ spp."
    ##### _Vibrio_ spp. {#vibrio}
    ??? task "`SRST2`: Vibrio characterization ==_for Illumina only_=="

        The `SRST2 Vibrio characterization` task detects sequences for *Vibrio* spp. characterization using Illumina sequence reads and a database of target sequence that are traditionally used in PCR methods. The sequences included in the database are as follows:
        
        | Sequence name | Sequence role | Purpose in database |
        | --- | --- | --- |
        | *toxR* | Transcriptional activator | Species marker where presence identifies *V. cholerae*  |
        | *ompW* | Outer Membrane Protein | Species marker where presence identifies *V. cholerae*  |
        | *ctxA* | Cholera toxin | Indicates cholera toxin production |
        | *tcpA*_classical | Toxin co-pilus A allele associated with the Classical biotype | Used to infer identity as Classical biotype |
        | tcpA_ElTor | Toxin co-pilus A allele associated with the El Tor biotype | Used to infer identity as El Tor biotype |
        | *wbeN* | O antigen encoding region | Used to infer identity as O1 serogroup |
        | *wbfR* | O antigen encoding region | Used to infer identity as O139 serogroup |

        !!! techdetails "SRST2 Technical Details"

            |  | Links |
            | --- | --- |
            | Task | [task_srst2_vibrio.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/vibrio/task_srst2_vibrio.wdl) |
            | Software Source Code | [srst2](https://github.com/katholt/srst2) |
            | Software Documentation | [srst2](https://github.com/katholt/srst2) |
            | Database Description | [Docker container](https://github.com/StaPH-B/docker-builds/tree/master/build-files/srst2/0.2.0-vibrio-230224) |
    
    ??? task "`Abricate`: Vibrio characterization"
        
        The `Abricate` Vibrio characterization task detects sequences for *Vibrio* spp. characterization using genome assemblies and the abricate "vibrio" database. The sequences included in the database are as follows:
        
        | Sequence name | Sequence role | Purpose in database |
        | --- | --- | --- |
        | *toxR* | Transcriptional activator | Species marker where presence identifies *V. cholerae*  |
        | *ompW* | Outer Membrane Protein | Species marker where presence identifies *V. cholerae*  |
        | *ctxA* | Cholera toxin | Indicates cholera toxin production |
        | *tcpA*_classical | Toxin co-pilus A allele associated with the Classical biotype | Used to infer identity as Classical biotype |
        | tcpA_ElTor | Toxin co-pilus A allele associated with the El Tor biotype | Used to infer identity as El Tor biotype |
        | *wbeN* | O antigen encoding region | Used to infer identity as O1 serogroup |
        | *wbfR* | O antigen encoding region | Used to infer identity as O139 serogroup |
        
        !!! techdetails "Abricate Technical Details"
            |  | Links |
            | --- | --- |
            | Task | [task_abricate_vibrio.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/vibrio/task_srst2_vibrio.wdl) |
            | Software Source Code | [abricate](https://github.com/tseemann/abricate) |
            | Software Documentation | [abricate](https://github.com/tseemann/abricate) |
            | Database Description | [Docker container](https://github.com/StaPH-B/docker-builds/tree/master/build-files/abricate/1.0.1-vibrio-cholera) |
    
    ??? task "`Vibecheck`: Vibrio cholerae classificaiton"
        
        The `Vibecheck` task classifies _V. cholerae_ sequences into canonical lineages (T1-T17) using variant frequency demixing. 
        
        !!! techdetails "Vibecheck Technical Details"
            |  | Links |
            | --- | --- |
            | Task | [task_vibecheck_vibrio.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/vibrio/task_vibecheck_vibrio.wdl) |
            | Software Source Code | [Vibecheck](https://github.com/CholGen/Vibecheck) |
            | Software Documentation | [Vibecheck](https://github.com/CholGen/Vibecheck) |
            | Database Description | [Docker container](https://hub.docker.com/r/watronfire/vibecheck) |

### Outputs

<div class="searchable-table" markdown="1">

| **Variable** | **Type** | **Description** | **Workflow** |
|---|---|---|---|
| abricate_abaum_database | String | Database of reference A. baumannii plasmid typing genes used for plasmid typing | FASTA, ONT, PE, SE |
| abricate_abaum_docker | String | Docker file used for running abricate | FASTA, ONT, PE, SE |
| abricate_abaum_plasmid_tsv | File | <https://github.com/tseemann/abricate#output> containing a row for each A. baumannii plasmid type gene found in the sample | FASTA, ONT, PE, SE |
| abricate_abaum_plasmid_type_genes | String | A. baumannii Plasmid typing genes found in the sample; from GENE column in <https://github.com/tseemann/abricate#output> | FASTA, ONT, PE, SE |
| abricate_abaum_version | String | Version of abricate used for A. baumannii plasmid typing | FASTA, ONT, PE, SE |
| abricate_database | String | Database of reference used with Abricate | FASTA, ONT, PE, SE |
| abricate_docker | String | Docker file used for running abricate | FASTA, ONT, PE, SE |
| abricate_genes | String | Genes found in the sample; from GENE column in <https://github.com/tseemann/abricate#output> | FASTA, ONT, PE, SE |
| abricate_results_tsv | File | <https://github.com/tseemann/abricate#output> containing a row for each gene found in the sample | FASTA, ONT, PE, SE |
| abricate_version | String | Version of abricate used for A. baumannii plasmid typing | FASTA, ONT, PE, SE |
| abricate_vibrio_biotype | String | Biotype classification according to tcpA gene sequence (Classical or ElTor) | FASTA, ONT, PE, SE |
| abricate_vibrio_ctxA | String | Presence or absence of the ctxA gene | FASTA, ONT, PE, SE |
| abricate_vibrio_detailed_tsv | File | Detailed ABRicate output file | FASTA, ONT, PE, SE |
| abricate_vibrio_ompW | String | Presence or absence of the ompW gene | FASTA, ONT, PE, SE |
| abricate_vibrio_serogroup | String | Serotype classification as O1 (wbeN gene), O139 (wbfR gene) or not detected. | FASTA, ONT, PE, SE |
| abricate_vibrio_toxR | String | Presence or absence of the toxR gene | FASTA, ONT, PE, SE |
| abricate_vibrio_version | String | The abricate version run | FASTA, ONT, PE, SE |
| agrvate_agr_canonical | String | Canonical or non-canonical agrD | FASTA, ONT, PE, SE |
| agrvate_agr_group | String | Agr group | FASTA, ONT, PE, SE |
| agrvate_agr_match_score | String | Match score for agr group | FASTA, ONT, PE, SE |
| agrvate_agr_multiple | String | If multiple agr groups were found | FASTA, ONT, PE, SE |
| agrvate_agr_num_frameshifts | String | Number of frameshifts found in CDS of extracted agr operon | FASTA, ONT, PE, SE |
| agrvate_docker | String | The docker used for AgrVATE | FASTA, ONT, PE, SE |
| agrvate_results | File | A gzipped tarball of all results | FASTA, ONT, PE, SE |
| agrvate_summary | File | The summary file produced | FASTA, ONT, PE, SE |
| agrvate_version | String | The version of AgrVATE used | FASTA, ONT, PE, SE |
| amr_results_csv | File | CSV formatted AMR profile | FASTA, ONT, PE, SE |
| amr_results_pdf | File | PDF formatted AMR profile |  FASTA, ONT, PE, SE |
| amr_search_results | File | JSON formatted AMR profile including BLAST results | FASTA, ONT, PE, SE |
| amr_search_docker | String | Docker image used to run AMR_Search | FASTA, ONT, PE, SE |
| amr_search_version | String | Version of AMR_Search libraries used | FASTA, ONT, PE, SE |
| amrfinderplus_all_report | File | Output TSV file from AMRFinderPlus (described <https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus#fields>) | FASTA, ONT, PE, SE |
| amrfinderplus_amr_betalactam_betalactam_genes | String | Beta-lactam AMR genes identified by AMRFinderPlus that are known to confer resistance to beta-lactams  | FASTA, ONT, PE, SE |
| amrfinderplus_amr_betalactam_carbapenem_genes | String | Beta-lactam AMR genes identified by AMRFinderPlus that are known to confer resistance to carbapenem  | FASTA, ONT, PE, SE |
| amrfinderplus_amr_betalactam_cephalosporin_genes | String | Beta-lactam AMR genes identified by AMRFinderPlus that are known to confer resistance to cephalosporin | FASTA, ONT, PE, SE |
| amrfinderplus_amr_betalactam_cephalothin_genes | String | Beta-lactam AMR genes identified by AMRFinderPlus that are known to confer resistance to cephalothin | FASTA, ONT, PE, SE |
| amrfinderplus_amr_betalactam_genes | String | Beta-lactam AMR genes identified by AMRFinderPlus | FASTA, ONT, PE, SE |
| amrfinderplus_amr_betalactam_methicillin_genes | String | Beta-lactam AMR genes identified by AMRFinderPlus that are known to confer resistance to methicilin | FASTA, ONT, PE, SE |
| amrfinderplus_amr_classes | String | AMRFinderPlus predictions for classes of drugs that genes found in the reads are known to confer resistance to | FASTA, ONT, PE, SE |
| amrfinderplus_amr_core_genes | String | AMR genes identified by AMRFinderPlus where the scope is "core" | FASTA, ONT, PE, SE |
| amrfinderplus_amr_plus_genes | String | AMR genes identified by AMRFinderPlus where the scope is "plus" | FASTA, ONT, PE, SE |
| amrfinderplus_amr_report | File | TSV file detailing AMR genes only, from the amrfinderplus_all_report | FASTA, ONT, PE, SE |
| amrfinderplus_amr_subclasses | String | More specificity about the drugs that genes identified in the reads confer resistance to | FASTA, ONT, PE, SE |
| amrfinderplus_db_version | String | AMRFinderPlus database version used | FASTA, ONT, PE, SE |
| amrfinderplus_stress_genes | String | Stress genes identified by AMRFinderPlus | FASTA, ONT, PE, SE |
| amrfinderplus_stress_report | File | TSV file detailing stress genes only, from the amrfinderplus_all_report | FASTA, ONT, PE, SE |
| amrfinderplus_version | String | AMRFinderPlus version used | FASTA, ONT, PE, SE |
| amrfinderplus_virulence_genes | String | Virulence genes identified by AMRFinderPlus | FASTA, ONT, PE, SE |
| amrfinderplus_virulence_report | File | TSV file detailing virulence genes only, from the amrfinderplus_all_report | FASTA, ONT, PE, SE |
| ani_highest_percent | Float | Highest ANI between query and any given reference genome (top species match) | FASTA, ONT, PE, SE |
| ani_highest_percent_bases_aligned | Float | Percentage of bases aligned between query genome and top species match | FASTA, ONT, PE, SE |
| ani_mummer_docker | String | Docker image used to run the ANI_mummer task | FASTA, ONT, PE, SE |
| ani_mummer_version | String | Version of MUMmer used | FASTA, ONT, PE, SE |
| ani_output_tsv | File | Full output TSV from ani-m | FASTA, ONT, PE, SE |
| ani_top_species_match | String | Species of genome with highest ANI to query FASTA | FASTA, ONT, PE, SE |
| assembly_fasta | File | Fasta file outputted from Flye or digger_denovo (Spades, SKEKA, or Megahit) | ONT, PE, SE |
| assembly_length | Int | Length of assembly (total contig length) as determined by QUAST | FASTA, ONT, PE, SE |
| assembler | String | Assembler used in digger_denovo subworkflow | PE, SE | 
| assembler_version | String | Version of the assembler used | PE, SE |
| bakta_gbff | File | Genomic GenBank format annotation file | FASTA, ONT, PE, SE |
| bakta_gff3 | File | Generic Feature Format Version 3 file | FASTA, ONT, PE, SE |
| bakta_plot | File | Bakta plot output PNG file summarizing annotated genome features such as coding sequences, RNA genes, and hypothetical proteins. | FASTA, ONT, PE, SE |
| bakta_summary | File | Bakta summary output TXT file | FASTA, ONT, PE, SE |
| bakta_tsv | File | Annotations as simple human readable TSV | FASTA, ONT, PE, SE |
| bakta_version | String | Bakta version used | FASTA, ONT, PE, SE |
| bandage_plot | File | Image file (PNG) visualizing the Flye assembly graph generated by Bandage | ONT |
| bandage_version | String | Version of Bandage used | ONT |
| bbduk_docker | String | BBDuk docker image used  | PE, SE |
| busco_database | String | BUSCO database used | FASTA, ONT, PE, SE |
| busco_docker | String | BUSCO docker image used | FASTA, ONT, PE, SE |
| busco_report | File | A plain text summary of the results in BUSCO notation | FASTA, ONT, PE, SE |
| busco_results | String | BUSCO results (see relevant toggle in [this block](https://theiagen.github.io/public_health_bioinformatics/latest/workflows/genomic_characterization/theiaprok/#post-assembly-tasks-performed-for-all-taxa)) | FASTA, ONT, PE, SE |
| busco_version | String | BUSCO software version used | FASTA, ONT, PE, SE |
| bwa_version | String | Version of BWA software used | ONT |
| cg_pipeline_docker | String | Docker file used for running CG-Pipeline on cleaned reads | PE, SE |
| cg_pipeline_report_clean | File | TSV file of read metrics from clean reads, including average read length, number of reads, and estimated genome coverage | PE, SE |
| cg_pipeline_report_raw | File | TSV file of read metrics from raw reads, including average read length, number of reads, and estimated genome coverage | PE, SE |
| clockwork_decontaminated_read1 | File | Decontaminated forward reads by Clockwork | PE |
| clockwork_decontaminated_read2 | File | Decontaminated reverse reads by Clockwork | PE |
| combined_mean_q_clean | Float | Mean quality score for the combined clean reads | PE |
| combined_mean_q_raw | Float | Mean quality score for the combined raw reads | PE |
| combined_mean_readlength_clean | Float | Mean read length for the combined clean reads | PE |
| combined_mean_readlength_raw | Float | Mean read length for the combined raw reads | PE |
| contigs_fastg | File | Assembly graph if megahit used for genome assembly | PE |
| contigs_gfa | File | Assembly graph output generated by SPAdes (Illumina: PE, SE) or Flye (ONT), used to visualize and evaluate genome assembly results. | ONT, PE, SE |
| contigs_lastgraph | File | Assembly graph if velvet used for genome assembly | PE |
| dnaapler_version | String | Version of dnaapler used | ONT |
| ectyper_predicted_serotype | String | Serotype predicted by ECTyper | FASTA, ONT, PE, SE |
| ectyper_results | File | TSV file of evidence for ECTyper predicted serotype (see <https://github.com/phac-nml/ecoli_serotyping#report-format>) | FASTA, ONT, PE, SE |
| ectyper_version | String | Version of ECTyper used | FASTA, ONT, PE, SE |
| emmtypingtool_docker | String | Docker image for emm-typing-tool | PE |
| emmtypingtool_emm_type | String | emm-type predicted | PE |
| emmtypingtool_results_xml | File | XML file with emm-typing-tool resuls | PE |
| emmtypingtool_version | String | Version of emm-typing-tool used | PE |
| est_coverage_clean | Float | Estimated coverage calculated from clean reads and genome length | ONT, PE, SE |
| est_coverage_raw | Float | Estimated coverage calculated from raw reads and genome length | ONT, PE, SE |
| fastp_html_report | File | The HTML report made with fastp | PE, SE |
| fastp_version | String | Version of fastp software used | PE, SE |
| fastq_scan_clean1_json | File | JSON file output from `fastq-scan` containing summary stats about clean forward read quality and length | PE, SE |
| fastq_scan_clean2_json | File | JSON file output from `fastq-scan` containing summary stats about clean reverse read quality and length | PE |
| fastq_scan_num_reads_clean_pairs | String | Number of read pairs after cleaning as calculated by fastq_scan | PE |
| fastq_scan_num_reads_clean1 | Int | Number of forward reads after cleaning as calculated by fastq_scan | PE, SE |
| fastq_scan_num_reads_clean2 | Int | Number of reverse reads after cleaning as calculated by fastq_scan | PE |
| fastq_scan_num_reads_raw_pairs | String | Number of input read pairs calculated by fastq_scan | PE |
| fastq_scan_num_reads_raw1 | Int | Number of input forward reads calculated by fastq_scan | PE, SE |
| fastq_scan_num_reads_raw2 | Int | Number of input reverse reads calculated by fastq_scan | PE |
| fastq_scan_raw1_json | File | JSON file output from `fastq-scan` containing summary stats about raw forward read quality and length | PE, SE |
| fastq_scan_raw2_json | File | JSON file output from `fastq-scan` containing summary stats about raw reverse read quality and length | PE |
| fastq_scan_version | String | Version of fastq-scan software used | PE, SE |
| fastqc_clean1_html | File | Graphical visualization of clean forward read quality from fastqc to open in an internet browser | PE, SE |
| fastqc_clean2_html | File | Graphical visualization of clean reverse read quality from fastqc to open in an internet browser | PE |
| fastqc_docker | String | Docker container used with fastqc | PE, SE |
| fastqc_num_reads_clean_pairs | String | Number of read pairs after cleaning by fastqc | PE |
| fastqc_num_reads_clean1 | Int | Number of forward reads after cleaning by fastqc | PE, SE |
| fastqc_num_reads_clean2 | Int | Number of reverse reads after cleaning by fastqc | PE |
| fastqc_num_reads_raw_pairs | String | Number of input read pairs by fastqc | PE |
| fastqc_num_reads_raw1 | Int | Number of input reverse reads by fastqc | PE, SE |
| fastqc_num_reads_raw2 | Int | Number of input reverse reads by fastqc | PE |
| fastqc_raw1_html | File | Graphical visualization of raw forward read quality from fastqc to open in an internet browser | PE, SE |
| fastqc_raw2_html | File | Graphical visualization of raw reverse read qualityfrom fastqc to open in an internet browser | PE |
| fastqc_version | String | Version of fastqc software used | PE, SE |
| flye_version | String | Version of Flye software used | ONT |
| gambit_closest_genomes | File | CSV file listing genomes in the GAMBIT database that are most similar to the query assembly | FASTA, ONT, PE, SE |
| gambit_db_version | String | Version of GAMBIT used | FASTA, ONT, PE, SE |
| gambit_docker | String | GAMBIT docker file used | FASTA, ONT, PE, SE |
| gambit_predicted_taxon | String | Taxon predicted by GAMBIT | FASTA, ONT, PE, SE |
| gambit_predicted_taxon_rank | String | Taxon rank of GAMBIT taxon prediction | FASTA, ONT, PE, SE |
| gambit_report | File | GAMBIT report in a machine-readable format | FASTA, ONT, PE, SE |
| gambit_version | String | Version of GAMBIT software used | FASTA, ONT, PE, SE |
| genotyphi_final_genotype | String | Final genotype call from GenoTyphi | ONT, PE, SE |
| genotyphi_genotype_confidence | String | Confidence in the final genotype call made by GenoTyphi | ONT, PE, SE |
| genotyphi_mykrobe_json | File | JSON file of GenoTyphi output, described <https://github.com/katholt/genotyphi#explanation-of-columns-in-the-output> | ONT, PE, SE |
| genotyphi_report_tsv | File | TSV file of GenoTyphi output, described <https://github.com/katholt/genotyphi#explanation-of-columns-in-the-output> | ONT, PE, SE |
| genotyphi_species | String | Species call from Mykrobe, used to run GenoTyphi | ONT, PE, SE |
| genotyphi_st_probes_percent_coverage | Float | Percentage coverage to the Typhi MLST probes | ONT, PE, SE |
| genotyphi_version | String | Version of GenoTyphi used | ONT, PE, SE |
| hicap_docker | String | Docker image used for hicap | ONT, PE, SE |
| hicap_genes | String | cap genes identified. genes on different contigs delimited by;. truncation shown by trailing * | ONT, PE, SE |
| hicap_results_tsv | File | TSV file of hicap output | ONT, PE, SE |
| hicap_serotype | String | hicap serotype | ONT, PE, SE |
| hicap_version | String | hicap version used | ONT, PE, SE |
| kaptive_k_locus | String | Best matching K locus identified by Kaptive | FASTA, ONT, PE, SE |
| kaptive_k_type | String | Best matching K type identified by Kaptive | FASTA, ONT, PE, SE |
| kaptive_kl_confidence | String | Kaptive’s confidence in the KL match (see <https://github.com/katholt/Kaptive/wiki/Interpreting-the-results>) | FASTA, ONT, PE, SE |
| kaptive_oc_locus | String | Best matching K locus identified by Kaptive | FASTA, ONT, PE, SE |
| kaptive_ocl_confidence | String | Kaptive’s confidence in the OCL match (see <https://github.com/katholt/Kaptive/wiki/Interpreting-the-results>) | FASTA, ONT, PE, SE |
| kaptive_output_file_k | File | TSV <https://github.com/katholt/Kaptive/wiki/How-to-run#output-filesfrom> the K locus from Kaptive | FASTA, ONT, PE, SE |
| kaptive_output_file_oc | File | TSV <https://github.com/katholt/Kaptive/wiki/How-to-run#output-filesfrom> the OC locus from Kaptive | FASTA, ONT, PE, SE |
| kaptive_version | String | Version of Kaptive used | FASTA, ONT, PE, SE |
| kleborate_docker | String | Kleborate docker image used | FASTA, ONT, PE, SE |
| kleborate_genomic_resistance_mutations | String | Genomic resistance mutations identifies by Kleborate | FASTA, ONT, PE, SE |
| kleborate_key_resistance_genes | String | Key resistance genes identified by Kleborate | FASTA, ONT, PE, SE |
| kleborate_klocus | String | Best matching K locus identified by  Kleborate via Kaptive | FASTA, ONT, PE, SE |
| kleborate_klocus_confidence | String | Kaptive’s confidence in the KL match (see <https://github.com/katholt/Kaptive/wiki/Interpreting-the-results>) | FASTA, ONT, PE, SE |
| kleborate_ktype | String | Best matching K type identified by  Kleborate via Kaptive | FASTA, ONT, PE, SE |
| kleborate_mlst_sequence_type | String | <https://github.com/katholt/Kleborate/wiki/MLST#multi-locus-sequence-typing-mlst> call by Kleborate | FASTA, ONT, PE, SE |
| kleborate_olocus | String | Best matching OC locus identified by  Kleborate via Kaptive | FASTA, ONT, PE, SE |
| kleborate_olocus_confidence | String | Kaptive’s confidence in the KL match (see <https://github.com/katholt/Kaptive/wiki/Interpreting-the-results>) | FASTA, ONT, PE, SE |
| kleborate_otype | String | Best matching OC type identified by  Kleborate via Kaptive | FASTA, ONT, PE, SE |
| kleborate_output_file | File | <https://github.com/katholt/Kleborate/wiki/Scores-and-counts> | FASTA, ONT, PE, SE |
| kleborate_resistance_score | String | Resistance score as given by kleborate | FASTA, ONT, PE, SE |
| kleborate_version | String | Version of Kleborate used | FASTA, ONT, PE, SE |
| kleborate_virulence_score | String | Virulence score as given by kleborate | FASTA, ONT, PE, SE |
| kmerfinder_database | String | Database used to run KmerFinder | FASTA, ONT, PE, SE |
| kmerfinder_docker | String | Docker image used to run KmerFinder | FASTA, ONT, PE, SE |
| kmerfinder_query_coverage | String | KmerFinder’s query coverage of the top hit result | FASTA, ONT, PE, SE |
| kmerfinder_results_tsv | File | Output TSV file created by KmerFinder | FASTA, ONT, PE, SE |
| kmerfinder_template_coverage | String |  | FASTA, ONT, PE, SE |
| kmerfinder_top_hit | String | Top hit species of KmerFinder | FASTA, ONT, PE, SE |
| kraken2_database | String | Kraken2 database used for the taxonomic assignment | ONT, PE, SE |
| kraken2_docker | String | Docker container for Kraken2 | ONT, PE, SE |
| kraken2_report | File | Report, in text format, of Kraken2 results | ONT, PE, SE |
| kraken2_version | String | Kraken2 version | ONT, PE, SE |
| legsta_predicted_sbt | String | Sequence based type predicted by Legsta | FASTA, ONT, PE, SE |
| legsta_results | File | TSV file of legsta results (see <https://github.com/tseemann/legsta#output>) | FASTA, ONT, PE, SE |
| legsta_version | String | Version of legsta used | FASTA, ONT, PE, SE |
| lissero_results | File | TSV results file from LisSero (see <https://github.com/MDU-PHL/LisSero#example-output>) | FASTA, ONT, PE, SE |
| lissero_serotype | String | Serotype predicted by LisSero | FASTA, ONT, PE, SE |
| lissero_version | String | Version of LisSero used | FASTA, ONT, PE, SE |
| medaka_model | String | Model used by Medaka | ONT |
| medaka_version | String | Version of Medaka used | ONT |
| meningotype_BAST | String | BAST type | FASTA, ONT, PE, SE |
| meningotype_FetA | String | FetA type | FASTA, ONT, PE, SE |
| meningotype_fHbp | String | fHbp type | FASTA, ONT, PE, SE |
| meningotype_NadA | String | NBA type | FASTA, ONT, PE, SE |
| meningotype_NHBA | String | NHBA type | FASTA, ONT, PE, SE |
| meningotype_PorA | String | PorA type | FASTA, ONT, PE, SE |
| meningotype_PorB | String | PorB type | FASTA, ONT, PE, SE |
| meningotype_serogroup | String | Serogroup | FASTA, ONT, PE, SE |
| meningotype_tsv | File | Full result file | FASTA, ONT, PE, SE |
| meningotype_version | String | Version of meningotype used | FASTA, ONT, PE, SE |
| midas_docker | String | MIDAS docker image used | PE, SE |
| midas_primary_genus | String | Genus of most abundant species in reads | PE, SE |
| midas_report | File | TSV report of full MIDAS results | PE, SE |
| midas_secondary_genus | String | Genus of the next most abundant species after removing all species of the most abundant genus | PE, SE |
| midas_secondary_genus_abundance | String | Relative abundance of secondary genus | PE, SE |
| midas_secondary_genus_coverage | String | Absolute coverage of secondary genus | PE, SE |
| n50_value | Int | N50 of assembly calculated by QUAST | FASTA, ONT, PE, SE |
| nanoplot_docker | String | Docker image for nanoplot | ONT |
| nanoplot_html_clean | File | Clean read file | ONT |
| nanoplot_html_raw | File | Raw read file | ONT |
| nanoplot_num_reads_clean1 | Int | Number of clean reads | ONT |
| nanoplot_num_reads_raw1 | Int | Number of raw reads | ONT |
| nanoplot_r1_est_coverage_clean | Float | Estimated coverage on the clean reads by nanoplot | ONT |
| nanoplot_r1_est_coverage_raw | Float | Estimated coverage on the raw reads by nanoplot | ONT |
| nanoplot_r1_mean_q_clean | Float | Mean quality score of clean forward reads | ONT |
| nanoplot_r1_mean_q_raw | Float | Mean quality score of raw forward reads | ONT |
| nanoplot_r1_mean_readlength_clean | Float | Mean read length of clean forward reads | ONT |
| nanoplot_r1_mean_readlength_raw | Float | Mean read length of raw forward reads | ONT |
| nanoplot_r1_median_q_clean | Float | Median quality score of clean forward reads | ONT |
| nanoplot_r1_median_q_raw | Float | Median quality score of raw forward reads | ONT |
| nanoplot_r1_median_readlength_clean | Float | Median read length of clean forward reads | ONT |
| nanoplot_r1_median_readlength_raw | Float | Median read length of raw forward reads | ONT |
| nanoplot_r1_n50_clean | Float | N50 of clean forward reads | ONT |
| nanoplot_r1_n50_raw | Float | N50 of raw forward reads | ONT |
| nanoplot_r1_stdev_readlength_clean | Float | Standard deviation read length of clean forward reads | ONT |
| nanoplot_r1_stdev_readlength_raw | Float | Standard deviation read length of raw forward reads | ONT |
| nanoplot_tsv_clean | File | Output TSV file created by nanoplot | ONT |
| nanoplot_tsv_raw | File | Output TSV file created by nanoplot | ONT |
| nanoplot_version | String | Version of nanoplot used for analysis | ONT |
| nanoq_version | String | Version of nanoq used in analysis | ONT |
| ngmaster_ngmast_porB_allele | String | porB allele number | FASTA, ONT, PE, SE |
| ngmaster_ngmast_sequence_type | String | NG-MAST sequence type | FASTA, ONT, PE, SE |
| ngmaster_ngmast_tbpB_allele | String | tbpB allele number | FASTA, ONT, PE, SE |
| ngmaster_ngstar_23S_allele | String | 23S rRNA allele number | FASTA, ONT, PE, SE |
| ngmaster_ngstar_gyrA_allele | String | gyrA allele number | FASTA, ONT, PE, SE |
| ngmaster_ngstar_mtrR_allele | String | mtrR allele number | FASTA, ONT, PE, SE |
| ngmaster_ngstar_parC_allele | String | parC allele number | FASTA, ONT, PE, SE |
| ngmaster_ngstar_penA_allele | String | penA allele number | FASTA, ONT, PE, SE |
| ngmaster_ngstar_ponA_allele | String | ponA allele number | FASTA, ONT, PE, SE |
| ngmaster_ngstar_porB_allele | String | porB allele number | FASTA, ONT, PE, SE |
| ngmaster_ngstar_sequence_type | String | NG-STAR sequence type | FASTA, ONT, PE, SE |
| ngmaster_tsv | File | TSV file with NG-MAST/NG-STAR typing | FASTA, ONT, PE, SE |
| ngmaster_version | String | ngmaster version | FASTA, ONT, PE, SE |
| number_contigs | Int | Total number of contigs in assembly | FASTA, ONT, PE, SE |
| pasty_all_serogroups | File | TSV file with details of each serogroup from pasty (see <https://github.com/rpetit3/pasty#example-prefixdetailstsv>) | FASTA, ONT, PE, SE |
| pasty_blast_hits | File | TSV file of BLAST hits from pasty (see <https://github.com/rpetit3/pasty#example-prefixblastntsv>) | FASTA, ONT, PE, SE |
| pasty_comment | String |  | FASTA, ONT, PE, SE |
| pasty_docker | String | pasty docker image used | FASTA, ONT, PE, SE |
| pasty_serogroup | String | Serogroup predicted by pasty | FASTA, ONT, PE, SE |
| pasty_serogroup_coverage | Float | The breadth of coverage of the O-antigen by pasty | FASTA, ONT, PE, SE |
| pasty_serogroup_fragments | Int | Number of BLAST hits included in the prediction (fewer is better) | FASTA, ONT, PE, SE |
| pasty_summary_tsv | File | TSV summary file of pasty outputs (see <https://github.com/rpetit3/pasty#example-prefixtsv>) | FASTA, ONT, PE, SE |
| pasty_version | String | Version of pasty used | FASTA, ONT, PE, SE |
| pbptyper_docker | String | pbptyper docker image used | FASTA, ONT, PE, SE |
| pbptyper_pbptype_predicted_tsv | File | TSV file of pbptyper results (see <https://github.com/rpetit3/pbptyper#example-prefixtsv>) | FASTA, ONT, PE, SE |
| pbptyper_predicted_1A_2B_2X | String | PBP type predicted by pbptyper | FASTA, ONT, PE, SE |
| pbptyper_version | String | Version of pbptyper used | FASTA, ONT, PE, SE |
| plasmidfinder_db_version | String | Version of PlasmidFnder used | FASTA, ONT, PE, SE |
| plasmidfinder_docker | String | PlasmidFinder docker image used | FASTA, ONT, PE, SE |
| plasmidfinder_plasmids | String | Names of plasmids identified by PlasmidFinder | FASTA, ONT, PE, SE |
| plasmidfinder_results | File | Output file from PlasmidFinder in TSV format | FASTA, ONT, PE, SE |
| plasmidfinder_seqs | File | Hit_in_genome_seq.fsa file produced by PlasmidFinder | FASTA, ONT, PE, SE |
| polypolish_version | String | Version of Polypolish used | ONT|
| poppunk_docker | String | PopPUNK docker image with GPSC database used | FASTA, ONT, PE, SE |
| poppunk_gps_cluster | String | GPS cluster predicted by PopPUNK  | FASTA, ONT, PE, SE |
| poppunk_GPS_db_version | String | Version of GPSC database used | FASTA, ONT, PE, SE |
| poppunk_gps_external_cluster_csv | File | GPSC v6 scheme designations | FASTA, ONT, PE, SE |
| poppunk_version | String | Version of PopPUNK used | FASTA, ONT, PE, SE |
| porechop_version | String | Version of Porechop used | ONT |
| prokka_gbk | File | GenBank file produced from Prokka annotation of input FASTA | FASTA, ONT, PE, SE |
| prokka_gff | File | Prokka output GFF3 file containing sequence and annotation (you can view this in IGV) | FASTA, ONT, PE, SE |
| prokka_sqn | File | A Sequin file for GenBank submission | FASTA, ONT, PE, SE |
| qc_check | String | A string that indicates whether or not the sample passes a set of pre-determined and user-provided QC thresholds | FASTA, ONT, PE, SE |
| qc_standard | File | The user-provided file that contains the QC thresholds used for the QC check | FASTA, ONT, PE, SE |
| quast_gc_percent | Float | The GC percent of your sample | FASTA, ONT, PE, SE |
| quast_report | File | TSV report from QUAST | FASTA, ONT, PE, SE |
| quast_version | String | Software version of QUAST used | FASTA, ONT, PE, SE |
| r1_mean_q_clean | Float | Mean quality score of clean forward reads | PE, SE |
| r1_mean_q_raw | Float | Mean quality score of raw forward reads | PE, SE |
| r1_mean_readlength_clean | Float | Mean read length of clean forward reads | PE, SE |
| r1_mean_readlength_raw | Float | Mean read length of raw forward reads | PE, SE |
| r2_mean_q_clean | Float | Mean quality score of clean reverse reads | PE |
| r2_mean_q_raw | Float | Mean quality score of raw reverse reads | PE |
| r2_mean_readlength_clean | Float | Mean read length of clean reverse reads | PE |
| r2_mean_readlength_raw | Float | Mean read length of raw reverse reads | PE |
| racon_version | String | Version of Racon used | ONT |
| rasusa_version | String | Version of RASUSA used for analysis | ONT |
| read_screen_raw | String | PASS or FAIL result from raw read screening; FAIL accompanied by the reason(s) for failure | ONT, PE, SE |
| read_screen_raw_tsv | File | Raw read screening report TSV depicting read counts, total read base pairs, and estimated genome length | ONT, PE, SE |
| read_screen_clean | String | PASS or FAIL result from clean read screening; FAIL accompanied by the reason(s) for failure | ONT, PE, SE |
| read_screen_clean_tsv | File | Clean read screening report TSV depicting read counts, total read base pairs, and estimated genome length | ONT, PE, SE |
| read1_clean | File | Clean forward reads file | ONT, PE, SE |
| read2_clean | File | Clean reverse reads file | PE |
| resfinder_db_version | String | Version of ResFinder database | FASTA, ONT, PE, SE |
| resfinder_docker | String | ResFinder docker image used | FASTA, ONT, PE, SE |
| resfinder_pheno_table | File | Table containing al AMR phenotypes | FASTA, ONT, PE, SE |
| resfinder_pheno_table_species | File | Table with species-specific AMR phenotypes | FASTA, ONT, PE, SE |
| resfinder_pointfinder_pheno_table | File | TSV showing presence(1)/absence(0) of predicted resistance against an antibiotic class | FASTA, ONT, PE, SE |
| resfinder_pointfinder_results | File | Predicted point mutations, grouped by the gene they occur in | FASTA, ONT, PE, SE |
| resfinder_predicted_pheno_resistance | String | Semicolon delimited list of antimicrobial drugs and associated genes and/or point mutations. <drug1>: <gene1>, <gene1>, <point_mutation1>; <drug2>: <gene3>, <gene4>; | FASTA, ONT, PE, SE |
| resfinder_predicted_resistance_Amp | String | States either Resistance or No Resistance predicted to Ampicillin based on resfinder phenotypic predictions | FASTA, ONT, PE, SE |
| resfinder_predicted_resistance_Axo | String | States either Resistance or No Resistance predicted to Ceftriaxone based on resfinder phenotypic predictions | FASTA, ONT, PE, SE |
| resfinder_predicted_resistance_Azm | String | States either Resistance or No Resistance predicted to Azithromycin based on resfinder phenotypic predictions | FASTA, ONT, PE, SE |
| resfinder_predicted_resistance_Cip | String | States either Resistance or No Resistance predicted to Ciprofloxacin based on resfinder phenotypic predictions | FASTA, ONT, PE, SE |
| resfinder_predicted_resistance_Smx | String | States either Resistance or No Resistance predicted to Sulfamethoxazole based on resfinder phenotypic predictions | FASTA, ONT, PE, SE |
| resfinder_predicted_resistance_Tmp | String | States either Resistance or No Resistance predicted to Trimothoprim based on resfinder phenotypic predictions | FASTA, ONT, PE, SE |
| resfinder_predicted_xdr_shigella | String | Final prediction of XDR Shigella status based on CDC definition. Explanation can be found in the description above this table. | FASTA, ONT, PE, SE |
| resfinder_results | File | Predicted resistance genes grouped by antibiotic class | FASTA, ONT, PE, SE |
| resfinder_seqs | File | FASTA of resistance gene sequences from user’s input sequence | FASTA, ONT, PE, SE |
| seq_platform | String | Sequencing platform input by the user | FASTA, ONT, PE, SE |
| seqsero2_note | String | Additional notes produced by SeqSero2 | FASTA, ONT, PE, SE |
| seqsero2_predicted_antigenic_profile | String | Antigenic profile predicted for Salmonella spp. by SeqSero2 | FASTA, ONT, PE, SE |
| seqsero2_predicted_contamination | String | Indicates whether contamination between Salmonella with different serotypes was predicted by SeqSero2 | FASTA, ONT, PE, SE |
| seqsero2_predicted_serotype | String | Serotype predicted by SeqSero2 | FASTA, ONT, PE, SE |
| seqsero2_report | File | TSV report produced by SeqSero2 | FASTA, ONT, PE, SE |
| seqsero2_version | String | Version of SeqSero2 used | FASTA, ONT, PE, SE |
| seroba_ariba_identity | String | Percentage identity between the query sequence and ARIBA-predicted serotype | PE |
| seroba_ariba_serotype | String | Serotype predicted by ARIBA, via SeroBA | PE |
| seroba_details | File | Detailed TSV file from SeroBA | PE |
| seroba_docker | String | SeroBA docker image used | PE |
| seroba_serotype | String | Serotype predicted by SeroBA | PE |
| seroba_version | String | SeroBA version used | PE |
| serotypefinder_docker | String | SerotypeFinder docker image used | FASTA, ONT, PE, SE |
| serotypefinder_report | File | TSV report produced by SerotypeFinder | FASTA, ONT, PE, SE |
| serotypefinder_serotype | String | Serotype predicted by SerotypeFinder | FASTA, ONT, PE, SE |
| shigatyper_docker | String | ShigaTyper docker image used | ONT, PE, SE |
| shigatyper_hits_tsv | File | Detailed TSV report from ShigaTyper (see <https://github.com/CFSAN-Biostatistics/shigatyper#example-prefix-hitstsv>) | ONT, PE, SE |
| shigatyper_ipaB_presence_absence | String | Presence (+) or absence (-) of ipaB identified by ShigaTyper | ONT, PE, SE |
| shigatyper_notes | String | Any notes output from ShigaTyper | ONT, PE, SE |
| shigatyper_predicted_serotype | String | Serotype predicted by ShigaTyper | ONT, PE, SE |
| shigatyper_summary_tsv | File | TSV summary report from ShigaTyper (see <https://github.com/CFSAN-Biostatistics/shigatyper#example-prefixtsv>) | ONT, PE, SE |
| shigatyper_version | String | Version of ShigaTyper used | ONT, PE, SE |
| shigeifinder_cluster | String | Shigella/EIEC cluster identified by ShigEiFinder | FASTA, ONT, PE, SE |
| shigeifinder_cluster_reads | String | Shigella/EIEC cluster identified by ShigEiFinder using read files as inputs | PE, SE |
| shigeifinder_docker | String | ShigEiFinder docker image used | FASTA, ONT, PE, SE |
| shigeifinder_docker_reads | String | ShigEiFinder docker image used using read files as inputs | PE, SE |
| shigeifinder_H_antigen | String | H-antigen gene identified by ShigEiFinder | FASTA, ONT, PE, SE |
| shigeifinder_H_antigen_reads | String | H-antigen gene identified by ShigEiFinder using read files as inputs | PE, SE |
| shigeifinder_ipaH_presence_absence | String | Presence (+) or absence (-) of ipaH identified by ShigEiFinder | FASTA, ONT, PE, SE |
| shigeifinder_ipaH_presence_absence_reads | String | Presence (+) or absence (-) of ipaH identified by ShigEiFinder using read files as inputs | PE, SE |
| shigeifinder_notes | String | Any notes output from ShigEiFinder | FASTA, ONT, PE, SE |
| shigeifinder_notes_reads | String | Any notes output from ShigEiFinder using read files as inputs | PE, SE |
| shigeifinder_num_virulence_plasmid_genes | String | Number of virulence plasmid genes identified by ShigEiFinder | FASTA, ONT, PE, SE |
| shigeifinder_num_virulence_plasmid_genes_reads | String | Number of virulence plasmid genes identified by ShigEiFinder using read files as inputs | PE, SE |
| shigeifinder_O_antigen | String | O-antigen gene identified by ShigEiFinder | FASTA, ONT, PE, SE |
| shigeifinder_O_antigen_reads | String | O-antigen gene identified by ShigEiFinder using read files as inputs | PE, SE |
| shigeifinder_report | File | TSV report from ShigEiFinder (see <https://github.com/LanLab/ShigEiFinder#shigeifinder>) | FASTA, ONT, PE, SE |
| shigeifinder_report_reads | File | TSV report from ShigEiFinder (see <https://github.com/LanLab/ShigEiFinder#shigeifinder>) using read files as inputs | PE, SE |
| shigeifinder_serotype | String | Serotype predicted by ShigEiFinder | FASTA, ONT, PE, SE |
| shigeifinder_serotype_reads | String | Serotype predicted by ShigEiFinder using read files as inputs | PE, SE |
| shigeifinder_version | String | ShigEiFinder version used | FASTA, ONT, PE, SE |
| shigeifinder_version_reads | String | ShigEiFinder version used using read files as inputs | PE, SE |
| sistr_allele_fasta | File | FASTA file of novel cgMLST alleles from SISTR | FASTA, ONT, PE, SE |
| sistr_allele_json | File | JSON file of cgMLST allele sequences and information (see <https://github.com/phac-nml/sistr_cmd#cgmlst-allele-search-results>) | FASTA, ONT, PE, SE |
| sistr_antigenic_formula | String | A field that aggregates the O, H1, and H2, antigen values in a single location for convenience | FASTA, ONT, PE, SE |
| sistr_cgmlst | File | CSV file of the cgMLST allelic profile from SISTR (see <https://github.com/phac-nml/sistr_cmd#cgmlst-allelic-profiles-output---cgmlst-profiles-cgmlst-profilescsv>) | FASTA, ONT, PE, SE |
| sistr_h1_antigens | String | The predicted H1 antigen | FASTA, ONT, PE, SE |
| sistr_h2_antigens | String | The predicted H2 antigen | FASTA, ONT, PE, SE |
| sistr_o_antigens | String | The predicted O antigen | FASTA, ONT, PE, SE |
| sistr_predicted_serotype | String | Serotype predicted by SISTR | FASTA, ONT, PE, SE |
| sistr_results | File | TSV results file produced by SISTR (see <https://github.com/phac-nml/sistr_cmd#primary-results-output--o-sistr-results>) | FASTA, ONT, PE, SE |
| sistr_serogroup | String | Serogroup predicted by SISTR | FASTA, ONT, PE, SE |
| sistr_serotype_cgmlst | String | cgMLST of the serogroup prediicted by SISTR | FASTA, ONT, PE, SE |
| sistr_version | String | Version of SISTR used | FASTA, ONT, PE, SE |
| sonneityping_final_genotype | String | Final genotype call from Mykrobe, via sonneityper | ONT, PE, SE |
| sonneityping_final_report_tsv | File | Detailed TSV report from mykrobe, via sonneityper (see <https://github.com/katholt/sonneityping#example-output>) | ONT, PE, SE |
| sonneityping_genotype_confidence | String | Confidence in the final genotype call from sonneityper  | ONT, PE, SE |
| sonneityping_genotype_name | String | Human readable alias for genotype, where available provided by sonneityper | ONT, PE, SE |
| sonneityping_mykrobe_docker | String | sonneityping docker image used | ONT, PE, SE |
| sonneityping_mykrobe_report_csv | File | CSV report from mykrobe via sonneityper (see <https://github.com/Mykrobe-tools/mykrobe/wiki/AMR-prediction-output#csv-file>) | ONT, PE, SE |
| sonneityping_mykrobe_report_json | File | JSON report from mykrobe via sonneityper (see <https://github.com/Mykrobe-tools/mykrobe/wiki/AMR-prediction-output#json-file>) | ONT, PE, SE |
| sonneityping_mykrobe_version | String | Version of sonneityping used | ONT, PE, SE |
| sonneityping_species | String | Species call from Mykrobe via sonneityping | ONT, PE, SE |
| spatyper_docker | String | spatyper docker image used  | FASTA, ONT, PE, SE |
| spatyper_repeats | String | order of identified repeats | FASTA, ONT, PE, SE |
| spatyper_tsv | File | TSV report with spatyper results | FASTA, ONT, PE, SE |
| spatyper_type | String | spa type | FASTA, ONT, PE, SE |
| spatyper_version | String | spatyper version used | FASTA, ONT, PE, SE |
| srst2_vibrio_biotype | String | Biotype classification according to tcpA gene sequence (Classical or ElTor) | PE, SE |
| srst2_vibrio_ctxA | String | Presence or absence of the ctxA gene | PE, SE |
| srst2_vibrio_detailed_tsv | File | Detailed <https://github.com/katholt/srst2> output file | PE, SE |
| srst2_vibrio_ompW | String | Presence or absence of the ompW gene | PE, SE |
| srst2_vibrio_serogroup | String | Serotype classification as O1 (wbeN gene), O139 (wbfR gene) or not detected. | PE, SE |
| srst2_vibrio_toxR | String | Presence or absence of the toxR gene | PE, SE |
| srst2_vibrio_version | String | The SRST2 version run | PE, SE |
| staphopiasccmec_docker | String | staphopia-sccmec docker image used | FASTA, ONT, PE, SE |
| staphopiasccmec_hamming_distance_tsv | File | staphopia-sccmec version | FASTA, ONT, PE, SE |
| staphopiasccmec_results_tsv | File | sccmec types and mecA presence | FASTA, ONT, PE, SE |
| staphopiasccmec_types_and_mecA_presence | String | staphopia-sccmec Hamming distance file | FASTA, ONT, PE, SE |
| staphopiasccmec_version | String | staphopia-sccmec presence and absence TSV file | FASTA, ONT, PE, SE |
| stxtyper_all_hits | String | Comma-separated list of matches of all types. Includes complete, partial, frameshift, internal stop, and novel hits. List is de-duplicated so multiple identical hits are only listed once. For example if 5 partial stx2 hits are detected in the genome, only 1 "stx2" will be listed in this field. To view the potential subtype for each partial hit, the user will need to view the stxtyper_report TSV file. | FASTA, ONT, PE, SE |
| stxtyper_ambiguous_hits | String | Comma-separated list of matches that have the OPERON output of "AMBIGUOUS". Ambiguous bases found in the query sequence (e.g., N) | FASTA, ONT, PE, SE |
| stxtyper_complete_operons | String | Comma-separated list of all COMPLETE operons detected by StxTyper. Show multiple hits if present in results. | FASTA, ONT, PE, SE |
| stxtyper_docker | String | Name of docker image used by the stxtyper task. | FASTA, ONT, PE, SE |
| stxtyper_extended_operons | String | Comma-separated list of all EXTENDED operons detected by StxTyper if coding sequence extends beyond the reference stop codon for one or both of the reference proteins. | FASTA, ONT, PE, SE |
| stxtyper_novel_hits | String | Comma-separated list of matches that have the OPERON output of "COMPLETE_NOVEL". Possible outputs "stx1", "stx2", or "stx1,stx2" | FASTA, ONT, PE, SE |
| stxtyper_num_hits | Int | Number of "hits" or rows present in the `stxtyper_report` TSV file | FASTA, ONT, PE, SE |
| stxtyper_partial_hits | String | Possible outputs "stx1", "stx2", or "stx1,stx2". Tells the user that there was a partial hit to either the A or B subunit, but does not describe which subunit, only the possible types from the PARTIAL matches. | FASTA, ONT, PE, SE |
| stxtyper_report | File | Raw results TSV file produced by StxTyper | FASTA, ONT, PE, SE |
| stxtyper_stx_frameshifts_or_internal_stop_hits | String | Comma-separated list of matches that have the OPERON output of "FRAMESHIFT" or "INTERNAL_STOP". Possible outputs "stx1", "stx2", or "stx1,stx2" | FASTA, ONT, PE, SE |
| stxtyper_version | String | Version of StxTyper used | FASTA, ONT, PE, SE |
| taxon_table_status | String | Status of the taxon table upload | FASTA, ONT, PE, SE |
| tbp_parser_average_genome_depth | Float | Optional output. Average genome depth across the reference genome  | ONT, PE, SE |
| tbp_parser_coverage_report | File | Optional output. TSV file with breadth of coverage of each gene associated with antimicrobial resistance in mycobacterium tuberculosis.  | ONT, PE, SE |
| tbp_parser_docker | String | Optional output. The docker image for tbp-parser | ONT, PE |
| tbp_parser_genome_percent_coverage | Float | Optional output. The percent of the genome covered at a depth greater than the specified minimum (default 10) | ONT, PE, SE |
| tbp_parser_laboratorian_report_csv | File | Optional output. Human-readable laboratorian report file containing the list of mutations found to be conferring resistance, both by WHO classification and expert rule implementation. The file contains the following columns: sample_id, tbprofiler_gene_name, tbprofiler_variant_locus_tag, tbprofiler_variant_substitution_type, tbprofiler_variant_substitution_nt, tbprofiler_variant_substitution_aa, confidence according to WHO, antimicrobial, depth, frequency, read_support, rationale ( WHO or expert rule), and warning if the coverage is below specified minimum (default 10) | ONT, PE, SE |
| tbp_parser_lims_report_csv | File | Optional output. LIMS digestable CSV report containing information on resistance for a set of antimicrobials ( No resistance to X detected, The detected genetic determinant(s) have uncertain significance, resistance to X cannot be ruled out and Genetic determinant(s) associated with resistance to X detected). For each antimicrobial, the mutations found are reported in the mutation_nucleotide; (mutation_protein) format, otherwise No mutations detected is reported.   | ONT, PE, SE |
| tbp_parser_looker_report_csv | File | Optional output. Looker digestible CSV report containing information on resistance for a set of antimicrobials (R for resistant, S for susceptible) | ONT, PE, SE |
| tbp_parser_version | String | Optional output. The version of tbp-parser | ONT, PE |
| tbprofiler_dr_type | String | Drug resistance type predicted by TB-Profiler (sensitive, Pre-MDR, MDR, Pre-XDR, XDR) | ONT, PE, SE |
| tbprofiler_main_lineage | String | Lineage(s) predicted by TBProfiler | ONT, PE, SE |
| tbprofiler_median_depth | Int | The median depth of the H37Rv TB reference genome covered by the sample | ONT, PE |
| tbprofiler_output_bai | File | Index BAM file generated by mapping sequencing reads to reference genome by TBProfiler | ONT, PE, SE |
| tbprofiler_output_bam | File | BAM alignment file produced by TBProfiler | ONT, PE, SE |
| tbprofiler_output_file | File | CSV report from TBProfiler | ONT, PE, SE |
| tbprofiler_output_vcf | File | VCF file output from TBProfiler; the concatenation of all of the different VCF files produced during TBProfiler analysis | ONT, PE, SE |
| tbprofiler_pct_reads_mapped | Float | The percentage of reads mapped to the H37Rv TB reference genome | ONT, PE |
| tbprofiler_resistance_genes | String | List of resistance mutations detected by TBProfiler | ONT, PE, SE |
| tbprofiler_sub_lineage | String | Sub-lineage(s) predicted by TBProfiler | ONT, PE, SE |
| tbprofiler_version | String | Version of TBProfiler used | ONT, PE, SE |
| theiaprok_fasta_analysis_date | String | Date of TheiaProk FASTA workflow execution | FASTA |
| theiaprok_fasta_version | String | Version of TheiaProk FASTA workflow execution | FASTA |
| theiaprok_illumina_pe_analysis_date | String | Date of TheiaProk PE workflow execution | PE |
| theiaprok_illumina_pe_version | String | Version of TheiaProk PE workflow execution | PE |
| theiaprok_illumina_se_analysis_date | String | Date of TheiaProk SE workflow execution | SE |
| theiaprok_illumina_se_version | String | Version of TheiaProk SE workflow execution | SE |
| theiaprok_ont_analysis_date | String | Date of TheiaProk ONT workflow execution | ONT |
| theiaprok_ont_version | String | Version of TheiaProk ONT workflow execution | ONT |
| trimmomatic_docker | String | Docker image used for trimmomatic | PE, SE |
| trimmomatic_version | String | Version of trimmomatic used | PE, SE |
| ts_mlst_allelic_profile | String | Profile of MLST loci and allele numbers predicted by MLST | FASTA, ONT, PE, SE |
| ts_mlst_docker | String | Docker image used for MLST | FASTA, ONT, PE, SE |
| ts_mlst_novel_alleles | File | FASTA file containing nucleotide sequence of any alleles that are not in the MLST database used by TheiaProk | FASTA, ONT, PE, SE |
| ts_mlst_predicted_st | String | ST predicted by MLST | FASTA, ONT, PE, SE |
| ts_mlst_pubmlst_scheme | String | PubMLST scheme used byMLST | FASTA, ONT, PE, SE |
| ts_mlst_results | File | TSV report with detailed MLST profile, including <https://github.com/tseemann/mlst#missing-data> | FASTA, ONT, PE, SE |
| ts_mlst_version | String | Version of Torsten Seeman’s MLST tool used | FASTA, ONT, PE, SE |
| virulencefinder_docker | String | VirulenceFinder docker image used | FASTA, ONT, PE, SE |
| virulencefinder_hits | String | Virulence genes detected by VirulenceFinder | FASTA, ONT, PE, SE |
| virulencefinder_report_tsv | File | Output TSV file created by VirulenceFinder | FASTA, ONT, PE, SE |
| vibecheck_lineage_report | File | Output CSV file created by Vibecheck | PE |
| vibecheck_top_lineage | String | Most likely lineage assigned by Vibecheck | PE |
| vibecheck_confidence | Float | How confidence lineage assignment is. 0 - uncertain; 100 - Very certain | PE |
| vibecheck_classification_notes | String | Additional information provided by Vibecheck during classification | PE |
| vibecheck_version | String | Version of Vibecheck used | PE |
| vibecheck_docker | String | Docker image used by Vibecheck task | PE |

</div>
