# TheiaViral Workflow Series

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Genomic Characterization](../../workflows_overview/workflows_type.md/#genomic-characterization) | [Viral](../../workflows_overview/workflows_kingdom.md/#viral) | PHB vX.X.X | Yes | Sample-level |

## TheiaViral Workflows

The **TheiaViral** workflows are designed for the assembly, quality assessment, and characterization of **non-amplicon-based** viral genomes. The TheiaViral workflows are particularly well-suited for handling diverse or recombinant pathogens, such as rabies virus and norovirus, which often pose challenges for traditional reference-based assembly methods. This is because existing references often fail to  adequately capture the genetic diversity present within a given sample. To address this, TheiaViral incorporates an unbiased reference selection step based on Average Nucleotide Identity (ANI). This step compares a de novo assembly of the sample against a comprehensive database of over 270,000 viral genomes, enabling selection of the most closely related reference. Currently, TheiaViral workflows perform variant calling and generate a consensus sequence, but additional downstream analyses are not yet implemented.

**What are the main differences between the TheiaViral and TheiaCov workflows?**

<div class="grid cards" markdown>

-   :material-database: **TheiaViral Workflows**

    ---

    * For non-amplicon-based viral genomes.
    * Supports relatively diverse and recombinant pathogens.
    * Utilizes an ANI-based reference selection from a *de novo* assembled genome.
    

-   :material-database: **TheiaCov Workflows**

    ---

    * For tiled PCR - amplicon-based viral genomes.
    * Supports a limited number of [pathogens](../../workflows/genomic_characterization/theiacov.md/#supported-organisms).
    * Utilizes a preselected manually curated reference genome.

</div>


There are currently two TheiaViral workflows designed to accommodate different kinds of input data:

1. Illumina paired-end sequencing (**TheiaViral_Illumina_PE**)
2. ONT sequencing (**TheiaViral_ONT**)


### Inputs

!!! dna ""
    ??? toggle "TheiaViral_Illumina_PE Input Read Data"

        The TheiaViral_Illumina_PE workflow takes in Illumina paired-end read data. Read file names should end with `.fastq` or `.fq`, with the optional addition of `.gz`. When possible, Theiagen recommends zipping files with [gzip](https://www.gnu.org/software/gzip/) before Terra uploads to minimize data upload time.

        By default, the workflow anticipates 2 x 150bp reads (i.e. the input reads were generated using a 300-cycle sequencing kit). Modifications to the optional parameter for `trim_minlen` may be required to accommodate shorter read data, such as the 2 x 75bp reads generated using a 150-cycle sequencing kit.

    ??? toggle "TheiaViral_ONT Input Read Data"

        The TheiaCoV_ONT workflow takes in base-called ONT read data. Read file names should end with `.fastq` or `.fq`, with the optional addition of `.gz`. When possible, Theiagen recommends zipping files with [gzip](https://www.gnu.org/software/gzip/) before uploading to Terra to minimize data upload time.

        **The ONT sequencing kit and base-calling approach can produce substantial variability in the amount and quality of read data. Genome assemblies produced by the TheiaCoV_ONT workflow must be quality assessed before reporting results.**

<div class="searchable-table" markdown="1">

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** | **Workflow** |
|---|---|---|---|---|---|---|
| theiaviral_ont | **read1** | File | Base-called ONT read file in FASTQ file format (compression optional) | | Required | ONT |
| theiaviral_ont | **samplename** | String | Name of the sample being analyzed | | Required | ONT |
| theiaviral_ont | **taxon** | String | Taxon ID or organism name of interest | | Required | ONT |
| bcftools_consensus | **cpu** | Int | Number of CPUs allocated for the task | 2 | Optional | ONT |
| bcftools_consensus | **disk_size** | Int | Disk size allocated for the task (in GB) | 100 | Optional | ONT |
| bcftools_consensus | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/staphb/bcftools:1.20" | Optional | ONT |
| bcftools_consensus | **memory** | Int | Memory allocated for the task (in GB) | 4 | Optional | ONT |
| checkv_consensus | **checkv_db** | File | CheckV database file | "gs://theiagen-large-public-files-rp/terra/databases/checkv/checkv-db-v1.5.tar.gz" | Optional | ONT |
| checkv_consensus | **cpu** | Int | Number of CPUs allocated for the task | 2 | Optional | ONT |
| checkv_consensus | **disk_size** | Int | Disk size allocated for the task (in GB) | 100 | Optional | ONT |
| checkv_consensus | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/staphb/checkv:1.0.3" | Optional | ONT |
| checkv_consensus | **memory** | Int | Memory allocated for the task (in GB) | 8 | Optional | ONT |
| checkv_denovo | **checkv_db** | File | CheckV database file | "gs://theiagen-large-public-files-rp/terra/databases/checkv/checkv-db-v1.5.tar.gz" | Optional | ONT |
| checkv_denovo | **cpu** | Int | Number of CPUs allocated for the task | 2 | Optional | ONT |
| checkv_denovo | **disk_size** | Int | Disk size allocated for the task (in GB) | 100 | Optional | ONT |
| checkv_denovo | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/staphb/checkv:1.0.3" | Optional | ONT |
| checkv_denovo | **memory** | Int | Memory allocated for the task (in GB) | 8 | Optional | ONT |
| clair3 | **clair3_model** | String | Model to be used by Clair3 | "r941_prom_hac_g360+g422" | Optional | ONT |
| clair3 | **cpu** | Int | Number of CPUs allocated for the task | 4 | Optional | ONT |
| clair3 | **disable_phasing** | Boolean | True/False that determines if variants should be called without whatshap phasing in full alignment calling | TRUE | Optional | ONT |
| clair3 | **disk_size** | Int | Disk size allocated for the task (in GB) | 100 | Optional | ONT |
| clair3 | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/theiagen/clair3-extra-models:1.0.10" | Optional | ONT |
| clair3 | **enable_gvcf** | Boolean | True/False that determines if an additional GVCF output should generated | FALSE | Optional | ONT |
| clair3 | **enable_haploid_precise** | Boolean | True/False that determines haploid calling mode where only 1/1 is considered as a variant | TRUE | Optional | ONT |
| clair3 | **include_all_contigs** | Boolean | True/False that determines if all contigs should be included in the output | TRUE | Optional | ONT |
| clair3 | **indel_min_af** | Float | Minimum Indel AF required for a candidate variant | 0.08 | Optional | ONT |
| clair3 | **memory** | Int | Memory allocated for the task (in GB) | 8 | Optional | ONT |
| clair3 | **snp_min_af** | Float | Minimum SNP AF required for a candidate variant | 0.08 | Optional | ONT |
| clair3 | **variant_quality** | Int | If set, variants with >$qual will be marked PASS, or LowQual otherwise | 2 | Optional | ONT |
| clean_check_reads | **cpu** | Int | Number of CPUs allocated for the task | 1 | Optional | ONT |
| clean_check_reads | **disk_size** | Int | Disk size allocated for the task (in GB) | 100 | Optional | ONT |
| clean_check_reads | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/bactopia/gather_samples:2.0.2" | Optional | ONT |
| clean_check_reads | **memory** | Int | Memory allocated for the task (in GB) | 2 | Optional | ONT |
| consensus_qc | **cpu** | Int | Number of CPUs allocated for the task | 1 | Optional | ONT |
| consensus_qc | **disk_size** | Int | Disk size allocated for the task (in GB) | 100 | Optional | ONT |
| consensus_qc | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1" | Optional | ONT |
| consensus_qc | **memory** | Int | Memory allocated for the task (in GB) | 2 | Optional | ONT |
| fasta_utilities | **cpu** | Int | Number of CPUs allocated for the task | 1 | Optional | ONT |
| fasta_utilities | **disk_size** | Int | Disk size allocated for the task (in GB) | 10 | Optional | ONT |
| fasta_utilities | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/biocontainers/seqkit:2.4.0--h9ee0642_0" | Optional | ONT |
| fasta_utilities | **memory** | Int | Memory allocated for the task (in GB) | 2 | Optional | ONT |
| flye | **additional_parameters** | String | Additional parameters for Flye assembler | | Optional | ONT |
| flye | **asm_coverage** | Int | Reduced coverage for initial disjointig assembly | | Optional | ONT |
| flye | **cpu** | Int | Number of CPUs allocated for the task | 4 | Optional | ONT |
| flye | **disk_size** | Int | Disk size allocated for the task (in GB) | 100 | Optional | ONT |
| flye | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/staphb/flye:2.9.4" | Optional | ONT |
| flye | **flye_polishing_iterations** | Int | Number of polishing iterations | 1 | Optional | ONT |
| flye | **genome_length** | Int | Expected genome length for assembly - requires `asm_coverage` | | Optional | ONT |
| flye | **keep_haplotypes** | Boolean | True/False to prevent collapsing alternative haplotypes | FALSE | Optional | ONT |
| flye | **memory** | Int | Memory allocated for the task (in GB) | 32 | Optional | ONT |
| flye | **minimum_overlap** | Int | Minimum overlap between reads | | Optional | ONT |
| flye | **no_alt_contigs** | Boolean | True/False to disable alternative contig generation | FALSE | Optional | ONT |
| flye | **read_error_rate** | Float | Expected error rate in reads | | Optional | ONT |
| flye | **read_type** | String | Type of read data for Flye | "--nano-hq" | Optional | ONT |
| flye | **scaffold** | Boolean | True/False to enable scaffolding using graph | FALSE | Optional | ONT |
| mask_low_coverage | **cpu** | Int | Number of CPUs allocated for the task | 2 | Optional | ONT |
| mask_low_coverage | **disk_size** | Int | Disk size allocated for the task (in GB) | 100 | Optional | ONT |
| mask_low_coverage | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/staphb/bedtools:2.31.0" | Optional | ONT |
| mask_low_coverage | **memory** | Int | Memory allocated for the task (in GB) | 8 | Optional | ONT |
| metabuli | **cpu** | Int | Number of CPUs allocated for the task | 4 | Optional | ONT |
| metabuli | **disk_size** | Int | Disk size allocated for the task (in GB) | 100 | Optional | ONT |
| metabuli | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/theiagen/metabuli:1.1.0" | Optional | ONT |
| metabuli | **extract_unclassified** | Boolean | True/False that determines if unclassified reads should be extracted and combined with the taxon specific extracted reads | FALSE | Optional | ONT |
| metabuli | **memory** | Int | Memory allocated for the task (in GB) | 8 | Optional | ONT |
| metabuli | **metabuli_db** | File | Metabuli database file | "gs://theiagen-large-public-files-rp/terra/databases/metabuli/refseq_virus-v223.tar.gz" | Optional | ONT |
| metabuli | **min_cov** | Float | Minimum query coverage threshold (0.0 - 1.0) | 0.0 | Optional | ONT |
| metabuli | **min_score** | Float | Minimum sequenece similarity score (0.0 - 1.0) | 0.0 | Optional | ONT |
| metabuli | **min_sp_score** | Float | Minimum score for species- or lower-level classification | 0.0 | Optional | ONT |
| metabuli | **taxonomy_path** | File | Path to taxonomy file | "gs://theiagen-large-public-files-rp/terra/databases/metabuli/new_taxdump.tar.gz" | Optional | ONT |
| minimap2 | **cpu** | Int | Number of CPUs allocated for the task | 2 | Optional | ONT |
| minimap2 | **disk_size** | Int | Disk size allocated for the task (in GB) | 100 | Optional | ONT |
| minimap2 | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/staphb/minimap2:2.22" | Optional | ONT |
| minimap2 | **memory** | Int | Memory allocated for the task (in GB) | 8 | Optional | ONT |
| minimap2 | **query2** | File | Internal component. Do not modify | | Optional | ONT |
| nanoplot_clean | **cpu** | Int | Number of CPUs allocated for the task | 4 | Optional | ONT |
| nanoplot_clean | **disk_size** | Int | Disk size allocated for the task (in GB) | 100 | Optional | ONT |
| nanoplot_clean | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/staphb/nanoplot:1.40.0" | Optional | ONT |
| nanoplot_clean | **max_length** | Int | Maximum read length for plotting | 100000 | Optional | ONT |
| nanoplot_clean | **memory** | Int | Memory allocated for the task (in GB) | 16 | Optional | ONT |
| nanoplot_raw | **cpu** | Int | Number of CPUs allocated for the task | 4 | Optional | ONT |
| nanoplot_raw | **disk_size** | Int | Disk size allocated for the task (in GB) | 100 | Optional | ONT |
| nanoplot_raw | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/staphb/nanoplot:1.40.0" | Optional | ONT |
| nanoplot_raw | **max_length** | Int | Maximum read length for plotting | 100000 | Optional | ONT |
| nanoplot_raw | **memory** | Int | Memory allocated for the task (in GB) | 16 | Optional | ONT |
| nanoq | **cpu** | Int | Number of CPUs allocated for the task | 2 | Optional | ONT |
| nanoq | **disk_size** | Int | Disk size allocated for the task (in GB) | 100 | Optional | ONT |
| nanoq | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/biocontainers/nanoq:0.9.0--hec16e2b_1" | Optional | ONT |
| nanoq | **max_read_length** | Int | Maximum read length to keep | 100000 | Optional | ONT |
| nanoq | **max_read_qual** | Int | Maximum read quality to keep | 10 | Optional | ONT |
| nanoq | **memory** | Int | Memory allocated for the task (in GB) | 2 | Optional | ONT |
| nanoq | **min_read_length** | Int | Minimum read length to keep | 500 | Optional | ONT |
| nanoq | **min_read_qual** | Int | Minimum read quality to keep | 10 | Optional | ONT |
| ncbi_datasets | **cpu** | Int | Number of CPUs allocated for the task | 1 | Optional | ONT |
| ncbi_datasets | **disk_size** | Int | Disk size allocated for the task (in GB) | 50 | Optional | ONT |
| ncbi_datasets | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/staphb/ncbi-datasets:16.38.1" | Optional | ONT |
| ncbi_datasets | **include_gbff** | Boolean | True/False to include gbff files in the output | FALSE | Optional | ONT |
| ncbi_datasets | **include_gff3** | Boolean | True/False to include gff3 files in the output | FALSE | Optional | ONT |
| ncbi_datasets | **memory** | Int | Memory allocated for the task (in GB) | 4 | Optional | ONT |
| ncbi_identify | **cpu** | Int | Number of CPUs allocated for the task | 1 | Optional | ONT |
| ncbi_identify | **disk_size** | Int | Disk size allocated for the task (in GB) | 50 | Optional | ONT |
| ncbi_identify | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/staphb/ncbi-datasets:16.38.1" | Optional | ONT |
| ncbi_identify | **memory** | Int | Memory allocated for the task (in GB) | 4 | Optional | ONT |
| ncbi_scrub_se | **cpu** | Int | Number of CPUs allocated for the task | 4 | Optional | ONT |
| ncbi_scrub_se | **disk_size** | Int | Disk size allocated for the task (in GB) | 100 | Optional | ONT |
| ncbi_scrub_se | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/ncbi/sra-human-scrubber:2.2.1" | Optional | ONT |
| ncbi_scrub_se | **memory** | Int | Memory allocated for the task (in GB) | 8 | Optional | ONT |
| ncbi_taxon_summary | **cpu** | Int | Number of CPUs allocated for the task | 1 | Optional | ONT |
| ncbi_taxon_summary | **disk_size** | Int | Disk size allocated for the task (in GB) | 50 | Optional | ONT |
| ncbi_taxon_summary | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/staphb/ncbi-datasets:16.38.1" | Optional | ONT |
| ncbi_taxon_summary | **memory** | Int | Memory allocated for the task (in GB) | 4 | Optional | ONT |
| parse_mapping | **cpu** | Int | Number of CPUs allocated for the task | 2 | Optional | ONT |
| parse_mapping | **disk_size** | Int | Disk size allocated for the task (in GB) | 100 | Optional | ONT |
| parse_mapping | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/staphb/samtools:1.17" | Optional | ONT |
| parse_mapping | **memory** | Int | Memory allocated for the task (in GB) | 8 | Optional | ONT |
| porechop | **cpu** | Int | Number of CPUs allocated for the task | 4 | Optional | ONT |
| porechop | **disk_size** | Int | Disk size allocated for the task (in GB) | 100 | Optional | ONT |
| porechop | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/staphb/porechop:0.2.4" | Optional | ONT |
| porechop | **memory** | Int | Memory allocated for the task (in GB) | 16 | Optional | ONT |
| porechop | **trimopts** | String | Additional trimming options for Porechop | | Optional | ONT |
| quast_denovo | **cpu** | Int | Number of CPUs allocated for the task | 2 | Optional | ONT |
| quast_denovo | **disk_size** | Int | Disk size allocated for the task (in GB) | 100 | Optional | ONT |
| quast_denovo | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/staphb/quast:5.0.2" | Optional | ONT |
| quast_denovo | **memory** | Int | Memory allocated for the task (in GB) | 2 | Optional | ONT |
| quast_denovo | **min_contig_length** | Int | Lower threshold for a contig length in bp. Shorter contigs won’t be taken into account | 500 | Optional | ONT |
| rasusa | **bases** | String | Explicitly set the number of bases required e.g., 4.3kb, 7Tb, 9000, 4.1MB. If this option is given, --coverage and --genome-size are ignored | | Optional | ONT |
| rasusa | **coverage** | Float | The desired coverage to sub-sample the reads to. If --bases is not provided, this option and --genome-size are required | 250 | Optional | ONT |
| rasusa | **cpu** | Int | Number of CPUs allocated for the task | 4 | Optional | ONT |
| rasusa | **disk_size** | Int | Disk size allocated for the task (in GB) | 100 | Optional | ONT |
| rasusa | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/staphb/rasusa:2.1.0" | Optional | ONT |
| rasusa | **frac** | Float | Subsample to a fraction of the reads - e.g., 0.5 samples half the reads | | Optional | ONT |
| rasusa | **memory** | Int | Memory allocated for the task (in GB) | 8 | Optional | ONT |
| rasusa | **num** | Int | Subsample to a specific number of reads | | Optional | ONT |
| rasusa | **read2** | File | Second read file for paired-end data | | Optional | PE |
| rasusa | **seed** | Int | Random seed for reproducibility | | Optional | ONT |
| raven | **cpu** | Int | Number of CPUs allocated for the task | 4 | Optional | ONT |
| raven | **disk_size** | Int | Disk size allocated for the task (in GB) | 100 | Optional | ONT |
| raven | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/theiagen/raven:1.8.3" | Optional | ONT |
| raven | **memory** | Int | Memory allocated for the task (in GB) | 16 | Optional | ONT |
| raven | **raven_identity** | Float | Threshold for overlap between two reads in order to construct an edge between them | 0.0 | Optional | ONT |
| raven | **raven_opts** | Int | Additional parameters for Raven assembler | | Optional | ONT |
| raven | **raven_polishing_iterations** | Int | Number of polishing iterations | 2 | Optional | ONT |
| read_mapping_stats | **cpu** | Int | Number of CPUs allocated for the task | 2 | Optional | ONT |
| read_mapping_stats | **disk_size** | Int | Disk size allocated for the task (in GB) | 100 | Optional | ONT |
| read_mapping_stats | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/staphb/samtools:1.15" | Optional | ONT |
| read_mapping_stats | **memory** | Int | Memory allocated for the task (in GB) | 8 | Optional | ONT |
| skani | **cpu** | Int | Number of CPUs allocated for the task | 2 | Optional | ONT |
| skani | **disk_size** | Int | Disk size allocated for the task (in GB) | 100 | Optional | ONT |
| skani | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/staphb/skani:0.2.2" | Optional | ONT |
| skani | **memory** | Int | Memory allocated for the task (in GB) | 4 | Optional | ONT |
| skani | **skani_db** | File | Skani database file | "gs://theiagen-large-public-files-rp/terra/databases/skani/skani_db_20250325.tar" | Optional | ONT |
| theiaviral_ont | **call_porechop** | Boolean | True/False to trim adapters with porechop | FALSE | Optional | ONT |
| theiaviral_ont | **genome_length** | Int | Expected genome length of taxon of interest | | Optional | ONT |
| theiaviral_ont | **min_allele_freq** | Float | Minimum allele frequency required for a variant to populate the consensus sequence | 0.6 | Optional | ONT |
| theiaviral_ont | **min_depth** | Int | Minimum read depth required for a variant to populate the consensus sequence | 10 | Optional | ONT |
| theiaviral_ont | **min_map_quality** | Int | Minimum mapping quality required for read alignments | 20 | Optional | ONT |
| theiaviral_ont | **read_extraction_rank** | String | Taxonomic rank to use for read extraction | "family" | Optional | ONT |
| theiaviral_ont | **reference_fasta** | File | Reference genome in FASTA format | | Optional | ONT |
| theiaviral_ont | **skip_rasusa** | Boolean | True/False to skip read subsampling with Rasusa | FALSE | Optional | ONT |
| theiaviral_ont | **skip_raven** | Boolean | True/False to skip assembly with Raven and instead use Flye | FALSE | Optional | ONT |
| theiaviral_ont | **skip_screen** | Boolean | True/False to skip read screening check prior to analysis | FALSE | Optional | ONT |
| version_capture | **docker** | String | Docker image used for the task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional | ONT |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) | | Optional | ONT |


### Core Tasks

#### Illumina Data Core Tasks

#### ONT Data Core Tasks

### Outputs

!!! caption "TEMP PLACEHOLDER: TheiaViral_ONT Workflow Diagram"
    ![TheiaViral Workflow Diagram](../../assets/figures/TheiaViral_ONT.png)