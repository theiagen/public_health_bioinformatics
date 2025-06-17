# Clair3 Variants

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Phylogenetic Construction](../../workflows_overview/workflows_type.md/#phylogenetic-construction) | [Any taxa](../../workflows_overview/workflows_kingdom.md/#any-taxa) | PHB v3.0.0 | Yes | Sample-level |

## Clair3_Variants_ONT

The `Clair3_Variants` workflow processes Oxford Nanopore Technologies (ONT) sequencing data to identify genetic variations compared to a reference genome. It combines minimap2's long-read alignment capabilities with Clair3's deep learning-based variant calling, designed specifically for ONT data characteristics. The workflow first aligns raw reads to a reference genome using ONT-optimized parameters, processes these alignments into sorted and indexed BAM files, and then employs Clair3's specialized models to detect variants including single nucleotide polymorphisms (SNPs) and insertions/deletions (indels). If enabled, the workflow can also identify longer indels and generate genome-wide variant calls in gVCF format for downstream analysis.

!!! caption "Clair3_Variants Workflow Diagram"
    ![Clair3_Variants Workflow Diagram](../../assets/figures/Clair3_Variants_WF_Diagram.png)

!!! tip "Example Use Cases"
   - **Variant Discovery**: Identify genetic variations in ONT sequencing data compared to a reference genome
   - **SNP and Indel Detection**: Accurately detect both small variants and longer indels
   - **Population Studies**: Generate standardized variant calls suitable for population-level analyses


### Supported Clair3 Models {#supported-clair3-models}

| Model | Chemistry | Source |
|-------|-----------|---------|
| `r941_prom_sup_g5014` | R9.4.1 | Clair3 1.0.10 Release |
| `r941_prom_hac_g360+g422` | R9.4.1 | Clair3 1.0.10 Release |
| `r941_prom_hac_g238` | R9.4.1 | Clair3 1.0.10 Release |
| `r1041_e82_400bps_sup_v500` | R10.4.1 | [nanoporetech/rerio](https://github.com/nanoporetech/rerio?tab=readme-ov-file#clair3-models) |
| `r1041_e82_400bps_hac_v500` | R10.4.1 | [nanoporetech/rerio](https://github.com/nanoporetech/rerio?tab=readme-ov-file#clair3-models) |
| `r1041_e82_400bps_sup_v410` | R10.4.1 | [nanoporetech/rerio](https://github.com/nanoporetech/rerio?tab=readme-ov-file#clair3-models) |
| `r1041_e82_400bps_hac_v410` | R10.4.1 | [nanoporetech/rerio](https://github.com/nanoporetech/rerio?tab=readme-ov-file#clair3-models) |
| `ont` | Various | Legacy (Recommended for Guppy3 and Guppy4) |
| `ont_guppy2` | Various | Legacy (For Guppy2 data) |
| `ont_guppy5` | Various | Legacy (For Guppy5 data) |

!!! hint ""
    The latest models for ONT are downloaded from the [nanoporetech/rerio github](https://github.com/nanoporetech/rerio?tab=readme-ov-file#clair3-models). Please let us know if there is a model not included you would like to see added. 

### Inputs

!!! warning "Note on Haploid Settings"
    Several parameters are set by default for haploid genome analysis:

    - clair3_disable_phasing is set to `true` since phasing is not relevant for haploid genomes
    - clair3_include_all_contigs is set to `true` to ensure complete genome coverage
    - clair3_enable_haploid_precise is set to `true` to only consider homozygous variants (1/1), which is appropriate for haploid genomes

<div class="searchable-table" markdown="1">

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| clair3_variants_ont | **clair3_cpu** | Int | Number of CPUs to use | 4 | Optional |
| clair3_variants_ont | **clair3_disable_phasing** | Boolean | Disable whatshap phasing | true | Optional |
| clair3_variants_ont | **clair3_disk_size** | Int | Disk size in GB | 100 | Optional |
| clair3_variants_ont | **clair3_docker** | String | Docker container for task | us-docker.pkg.dev/general-theiagen/staphb/clair3:1.0.10 | Optional |
| clair3_variants_ont | **clair3_enable_gvcf** | Boolean | Output gVCF format | false | Optional |
| clair3_variants_ont | **clair3_enable_haploid_precise** | Boolean | Enable haploid precise calling, only 1/1 is considered as a variant | true | Optional |
| clair3_variants_ont | **clair3_enable_long_indel** | Boolean | Enable long indel calling | false | Optional |
| clair3_variants_ont | **clair3_include_all_contigs** | Boolean | Call variants on all contigs, should always be true for non-human samples | true | Optional |
| clair3_variants_ont | **clair3_memory** | Int | Memory allocation in GB | 8 | Optional |
| clair3_variants_ont | **clair3_model** | String | Model name for variant calling (see [supported models](#supported-clair3-models) for available options) | r941_prom_hac_g360+g422 | Optional |
| clair3_variants_ont | **clair3_variant_quality** | Int | Minimum variant quality score | 2 | Optional |
| clair3_variants_ont | **read1** | File | ONT sequencing reads in FASTQ format | | Required |
| clair3_variants_ont | **reference_genome_file** | File | Reference genome in FASTA format | | Required |
| clair3_variants_ont | **samplename** | String | Name of Samples | | Required |
| minimap2 | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| minimap2 | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| minimap2 | **docker** | String | Docker image used for this task. | "us-docker.pkg.dev/general-theiagen/staphb/minimap2:2.22" | Optional |
| minimap2 | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| minimap2 | **query2** | File | Internal component. Do not modify | None | Do not modify, Optional |
| sam_to_sorted_bam | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| sam_to_sorted_bam | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| sam_to_sorted_bam | **docker** | String | Docker image used for this task. | "us-docker.pkg.dev/general-theiagen/staphb/samtools:1.17" | Optional |
| sam_to_sorted_bam | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| samtools_faidx | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional |
| samtools_faidx | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| samtools_faidx | **docker** | String | Docker image used for this task. | "us-docker.pkg.dev/general-theiagen/staphb/samtools:1.17" | Optional |
| samtools_faidx | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| version_capture | **docker** | String | Docker container for versioning | us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0 | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional |

</div>

### Workflow Tasks

??? task "`minimap2`: Read Alignment"

    `minimap2` is used with long read specific parameters (-L --cs --MD flags) to align ONT reads to the reference genome. These specialized parameters are essential for proper handling of long read error profiles, generation of detailed alignment information, and improved mapping accuracy for long reads.

    !!! techdetails "minimap2 Technical Details"
      | | Links |
      |---|---|
      | Task | [task_minimap2.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/alignment/task_minimap2.wdl) |
      | Software Source Code | [minimap2 on GitHub](https://github.com/lh3/minimap2) |
      | Software Documentation | [minimap2](https://lh3.github.io/minimap2) |
      | Original Publication(s) | [Minimap2: pairwise alignment for nucleotide sequences](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778) |

??? task "`samtools`: BAM Processing"

    The bam processing step aligns files through several coordinate-based steps to prepare for variant calling. The task converts SAM format to BAM, sorts the BAM file by coordinate, and creates a BAM index file. This processed BAM is required for Clair3's variant calling pipeline.

    !!! techdetails "samtools Technical Details"
      | | Links |
      |---|---|
      | Task | [task_samtools.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/data_handling/task_parse_mapping.wdl) |
      | Software Source Code | [samtools on GitHub](https://github.com/samtools/samtools) |
      | Software Documentation | [samtools](https://www.htslib.org/doc/samtools.html) |
      | Original Publication(s) | [The Sequence Alignment/Map format and SAMtools](https://doi.org/10.1093/bioinformatics/btp352)<br>[Twelve Years of SAMtools and BCFtools](https://doi.org/10.1093/gigascience/giab008) |

??? task "`samtools faidx`: Reference Genome Indexing"

    `samtools faidx` creates necessary index files for the reference. This indexing step is    essential for enabling efficient random access to the reference sequence during variant calling.

    !!! techdetails "samtools Technical Details"
      | | Links |
      |---|---|
      | Task | [task_samtools.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/data_handling/task_parse_mapping.wdl) |
      | Software Source Code | [samtools on GitHub](https://github.com/samtools/samtools) |
      | Software Documentation | [samtools](https://www.htslib.org/doc/samtools.html) |
      | Original Publication(s) | [The Sequence Alignment/Map format and SAMtools](https://doi.org/10.1093/bioinformatics/btp352)<br>[Twelve Years of SAMtools and BCFtools](https://doi.org/10.1093/gigascience/giab008) |

??? task "`Clair3`: Variant Calling"

    `Clair3` performs deep learning-based variant detection using a multi-stage approach. The process begins with pileup-based calling for initial variant identification, followed by full-alignment analysis for comprehensive variant detection. Results are merged into a final high-confidence call set.

    The variant calling pipeline employs specialized neural networks trained on ONT data to accurately identify:
    - Single nucleotide variants (SNVs)
    - Small insertions and deletions (indels)
    - Structural variants

    !!! techdetails "Clair3 Technical Details"
      |  | Links |
      | --- | --- |
      | Task | [task_clair3.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_clair3.wdl) |
      | Software Source Code | [Clair3 on GitHub](https://github.com/HKU-BAL/Clair3) |
      | Software Documentation | [Clair3 Documentation](https://github.com/HKU-BAL/Clair3?tab=readme-ov-file#usage) |
      | Original Publication(s) | [Symphonizing pileup and full-alignment for deep learning-based long-read variant calling](https://doi.org/10.1101/2021.12.29.474431) |

### Outputs

<div class="searchable-table" markdown="1">

| **Variable** | **Type** | **Description** |
|---|---|---|
| aligned_bam | File | Sorted BAM file containing the minimap2 alignments of reads to the reference genome |
| aligned_bai | File | Index file for the aligned BAM |
| aligned_fai | File | Index file for the reference genome |
| clair3_docker_image | String | Version of the Docker container used for Clair3 variant calling |
| clair3_model_used | String | Name of the Clair3 model used for variant calling |
| clair3_variants_vcf | File | Final merged VCF file containing high-confidence variant calls, combining results from both pileup and full-alignment approaches |
| clair3_variants_gvcf | File | Optional genome VCF file containing information about all genomic positions, including non-variant sites |
| clair3_variants_wf_version | String | Version of the PHB workflow used |
| clair3_version | String | Clair3 Version being used |
| samtools_version | String | Version of samtools used for BAM processing |

</div>