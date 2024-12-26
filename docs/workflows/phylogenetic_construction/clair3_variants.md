# Clair3 Variants

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Phylogenetic Construction](../../workflows_overview/workflows_type.md/#phylogenetic-construction) | [Bacteria](../../workflows_overview/workflows_kingdom.md/#bacteria), [Mycotics](../../workflows_overview/workflows_kingdom.md#mycotics), [Viral](../../workflows_overview/workflows_kingdom.md/#viral) | PHB v2.3.0 | Yes | Sample-level |

## Clair3_Variants_ONT

The `Clair3_Variants` workflow processes Oxford Nanopore Technologies (ONT) sequencing data to identify genetic variations compared to a reference genome. It combines minimap2's long-read alignment capabilities with Clair3's deep learning-based variant calling, designed specifically for ONT data characteristics. The workflow first aligns raw reads to a reference genome using ONT-optimized parameters, processes these alignments into sorted and indexed BAM files, and then employs Clair3's specialized models to detect variants including single nucleotide polymorphisms (SNPs) and insertions/deletions (indels). If enabled, the workflow can also identify longer indels and generate genome-wide variant calls in gVCF format for downstream analysis.

!!! caption "Clair3_Variants Workflow Diagram"
   [Workflow diagram to be added]

!!! tip "Example Use Cases"
   - **Variant Discovery**: Identify genetic variations in ONT sequencing data compared to a reference genome
   - **SNP and Indel Detection**: Accurately detect both small variants and longer indels
   - **Population Studies**: Generate standardized variant calls suitable for population-level analyses

### Inputs

<div class="searchable-table" markdown="1">

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| clair3_variants_ont | **read1** | File | ONT sequencing reads in FASTQ format | | Required |
| clair3_variants_ont | **reference_genome_file** | File | Reference genome in FASTA format | | Required |
| clair3_variants_ont | **samplename** | String | Name of Samples | | Required |
| clair3_variants_ont | **clair3_model** | String | Model name for variant calling | r941_prom_hac_g360+g422 | Optional |
| clair3_variants | **variant_quality** | Int | Minimum variant quality score | 2 | Optional |
| clair3_variants | **enable_gvcf** | Boolean | Output gVCF format | false | Optional |
| clair3_variants | **enable_long_indel** | Boolean | Enable long indel calling | false | Optional |
| clair3_variants | **memory** | Int | Memory allocation in GB | 8 | Optional |
| clair3_variants | **cpu** | Int | Number of CPUs to use | 4 | Optional |
| clair3_variants | **disk_size** | Int | Disk size in GB | 100 | Optional |
| clair3_variants | **docker** | String | Docker container for task | us-docker.pkg.dev/general-theiagen/staphb/clair3:1.0.10 | Optional |
| version_capture | **docker** | String | Docker container for versioning | us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0 | Optional |

</div>

### Workflow Tasks

`Clair3_Variants` executes several key tasks in sequence to process ONT reads and identify variants:

1. **Minimap2 Alignment**: Uses minimap2 with long read specific parameters (-L --cs --MD flags) to align reads to the reference genome. These specialized parameters are essential for:
  - Proper handling of long read error profiles
  - Generation of detailed alignment information
  - Improved mapping accuracy for long reads

2. **SAM/BAM Processing**: Converts and processes alignment files through several steps:
  - Converts SAM format to BAM 
  - Sorts BAM file by coordinate
  - Creates BAM index file
  - This processed BAM is required for Clair3's variant calling

3. **Reference Genome Indexing**: Creates necessary index files for the reference genome using samtools faidx

4. **Clair3 Variant Calling**: Performs deep learning-based variant detection using:
  - Pileup-based calling for initial variant identification
  - Full-alignment analysis for comprehensive variant detection
  - Merges results into final high-confidence call set

!!! techdetails "Technical Details"
   |  | Links |
   | --- | --- |
   | Task | [task_clair3_variants.wdl](../tasks/gene_typing/variant_detection/task_clair3_variants.wdl) |
   | Software Source Code | [Clair3 on GitHub](https://github.com/HKU-BAL/Clair3) |
   | Software Documentation | [Clair3 Documentation](https://github.com/HKU-BAL/Clair3#clair3-documentation) |

### Outputs

<div class="searchable-table" markdown="1">

| **Variable** | **Type** | **Description** |
|---|---|---|
| aligned_bam | File | Sorted BAM file containing the minimap2 alignments of reads to the reference genome |
| aligned_bai | File | Index file for the aligned BAM |
| aligned_fai | File | Index file for the reference genome |
| clair3_variants_final_vcf | File | Final merged VCF file containing high-confidence variant calls, combining results from both pileup and full-alignment approaches |
| clair3_variants_pileup_vcf | File | VCF file containing variants identified using the pileup-based approach |
| clair3_variants_full_alignment_vcf | File | VCF file containing variants identified using the full-alignment analysis |
| clair3_variants_gvcf | File | Optional genome VCF file containing information about all genomic positions, including non-variant sites |
| clair3_docker_image | String | Version of the Docker container used for Clair3 variant calling |
| clair3_model_used | String | Name of the Clair3 model used for variant calling |
| samtools_version | String | Version of samtools used for BAM processing |
| clair3_variants_wf_version | String | Version of the Clair3 workflow used |

</div>

!!! warning "Note on VCF Outputs"
   The workflow produces multiple VCF files with different levels of filtering and analysis:
   - The `clair3_variants_final_vcf` represents the most confident calls and should be used for most analyses
   - The `pileup` and `full_alignment` VCFs can be useful for investigating specific variants or debugging
   - The gVCF output (if enabled) provides complete genomic context including non-variant positions