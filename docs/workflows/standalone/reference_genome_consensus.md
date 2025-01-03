# Reference_Genome_Consensus

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Standalone](../../workflows_overview/workflows_type.md/#standalone) | [Viral](../../workflows_overview/workflows_kingdom.md/#viral) | PHB vX.X.X | yes| Sample-level |

## Reference_Genome_Consensus_PHB

A workflow for reference-guided genome assembly, read trimming, alignment, coverage filtering, and consensus sequence generation from ONT basecalled FASTQ files.

### Inputs

<div class="searchable-table" markdown="1">

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| reference_genome_consensus | **fastq_file** | File | Input ONT basecalled FASTQ file | N/A | Required |
| reference_genome_consensus | **reference_file** | File | A file containing one or more reference genomes in FASTA format. Each genome must be properly formatted, with unique and descriptive headers. | N/A | Required |
| reference_genome_consensus | **sample_name** | String | Sample name for output file naming | N/A | Required |
| reference_genome_consensus | **left_trim** | Int | Number of bases to trim from the left | 25 | Optional |
| reference_genome_consensus | **right_trim** | Int | Number of bases to trim from the right | 25 | Optional |
| reference_genome_consensus | **min_length** | Int | Minimum read length to retain | 50 | Optional |
| reference_genome_consensus | **min_coverage** | Int | Minimum coverage for filtering regions | 10 | Optional |
| reference_genome_consensus | **medaka_model_override** | String | Override for Medaka model selection (optional) | N/A | Optional |
| fastqc_se | **cpu** | Int | Number of CPUs allocated for FastQC task | 2 | Optional |
| fastqc_se | **disk_size** | Int | Disk size allocated for FastQC task | 100 | Optional |
| fastqc_se | **docker** | String | Docker container for FastQC task | us-docker.pkg.dev/general-theiagen/staphb/fastqc:0.12.1 | Optional |
| fastqc_se | **memory** | Int | Memory allocated for FastQC task | 4 | Optional |
| seqtk_trim | **cpu** | Int | Number of CPUs allocated for Seqtk trimming task | 2 | Optional |
| seqtk_trim | **disk_size** | Int | Disk size allocated for Seqtk trimming task | 100 | Optional |
| seqtk_trim | **docker** | String | Docker container for Seqtk trimming task | us-docker.pkg.dev/general-theiagen/staphb/seqtk:1.4 | Optional |
| seqtk_trim | **memory** | Int | Memory allocated for Seqtk trimming task | 4 | Optional |
| plot_read_length_distribution | **cpu** | Int | Number of CPUs allocated for plotting read lengths | 1 | Optional |
| plot_read_length_distribution | **disk_size** | Int | Disk size allocated for plotting read lengths | 10 | Optional |
| plot_read_length_distribution | **docker** | String | Docker container for plotting read lengths | us-docker.pkg.dev/general-theiagen/theiagen/utility:1.2 | Optional |
| plot_read_length_distribution | **memory** | Int | Memory allocated for plotting read lengths | 4 | Optional |
| minimap2_align_ont | **cpu** | Int | Number of CPUs allocated for Minimap2 alignment | 2 | Optional |
| minimap2_align_ont | **disk_size** | Int | Disk size allocated for Minimap2 alignment | 100 | Optional |
| minimap2_align_ont | **docker** | String | Docker container for Minimap2 alignment | us-docker.pkg.dev/general-theiagen/staphb/minimap2:2.22 | Optional |
| minimap2_align_ont | **memory** | Int | Memory allocated for Minimap2 alignment | 8 | Optional |
| filter_coverage | **cpu** | Int | Number of CPUs allocated for filtering coverage | 4 | Optional |
| filter_coverage | **disk_size** | Int | Disk size allocated for filtering coverage | 50 | Optional |
| filter_coverage | **docker** | String | Docker container for filtering coverage | us-docker.pkg.dev/general-theiagen/staphb/bedtools:2.31.1 | Optional |
| filter_coverage | **memory** | Int | Memory allocated for filtering coverage | 8 | Optional |
| medaka_consensus | **auto_model** | Boolean | Automatically determine the Medaka model | true | Optional |
| medaka_consensus | **cpu** | Int | Number of CPUs allocated for Medaka consensus generation | 4 | Optional |
| medaka_consensus | **disk_size** | Int | Disk size allocated for Medaka consensus generation | 100 | Optional |
| medaka_consensus | **docker** | String | Docker container for Medaka consensus generation | us-docker.pkg.dev/general-theiagen/staphb/medaka:2.0.1 | Optional |
| medaka_consensus | **memory** | Int | Memory allocated for Medaka consensus generation | 16 | Optional |
| samtools_process | **cpu** | Int | Number of CPUs allocated for Samtools processing | 4 | Optional |
| samtools_process | **disk_size** | Int | Disk size allocated for Samtools processing | 50 | Optional |
| samtools_process | **docker** | String | Docker container for Samtools processing | us-docker.pkg.dev/general-theiagen/staphb/samtools:1.15 | Optional |
| samtools_process | **memory** | Int | Memory allocated for Samtools processing | 16 | Optional |

</div>

### Reference Genome File Format

The reference genome file `reference_file` must be in **FASTA format**. Each genome must be properly formatted, with unique and descriptive headers. For example:

```plaintext
>Genome_1
ATCGATCGATCG
>Genome_2
CGTAGCTAGCTA
```

The file must be in FASTA format (extensions such as .fa or .fasta).

Ensure the reference genomes are suitable for alignment with ONT reads and compatible with downstream tools like Minimap2.

### Workflow Tasks

This workflow includes the following tasks:

??? task "`fastqc`: Run FastQC"
    Quality control check of the input FASTQ file.

    !!! techdetails "FastQC Technical Details"
        |  | Links |
        | --- | --- |
        | Software Source Code | [FastQC](https://github.com/s-andrews/FastQC) |
        | Software Documentation | [FastQC Docs](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) |

??? task "`seqtk`: Trimming and Filtering Reads"
    Trims and filters reads based on user-defined parameters such as left trim, right trim, and minimum read length.

    !!! techdetails "Seqtk Technical Details"
        |  | Links |
        | --- | --- |
        | Software Source Code | [Seqtk](https://github.com/lh3/seqtk) |

??? task "`minimap2`: Align Reads to Reference"
    Aligns filtered reads to the reference genome.

    !!! techdetails "Minimap2 Technical Details"
        |  | Links |
        | --- | --- |
        | Software Source Code | [Minimap2](https://github.com/lh3/minimap2) |
        | Original Publication(s) | [Li 2018](https://doi.org/10.1093/bioinformatics/bty191) |

??? task "`samtools`: Process BAM File and Calculate Coverage"
    Converts SAM to BAM, sorts and indexes the BAM file, and calculates coverage.

    !!! techdetails "Samtools Technical Details"
        |  | Links |
        | --- | --- |
        | Software Source Code | [Samtools](http://www.htslib.org/) |

??? task "`filter_coverage`: Extract Regions with Sufficient Coverage"
    Filters aligned reads to retain regions meeting the minimum coverage threshold.

    !!! techdetails "Filter Coverage Technical Details"
        |  | Links |
        | --- | --- |
        | Software Source Code | [Bedtools](https://bedtools.readthedocs.io/) |

??? task "`medaka`: Generate Polished Consensus Sequence"
    Generates a consensus sequence using Medaka.

    !!! techdetails "Medaka Technical Details"
        |  | Links |
        | --- | --- |
        | Software Source Code | [Medaka](https://github.com/nanoporetech/medaka) |
        | Software Documentation | [Medaka Docs](https://nanoporetech.github.io/medaka/) |

### Outputs

| **Variable** | **Type** | **Description** |
|---|---|---|
| fastqc_report | File | HTML report from FastQC |
| read_length_plot | File | Plot of read length distribution |
| coverage_file | File | Coverage file generated by Samtools |
| medaka_consensus_fasta | File | Polished consensus sequence in FASTA format |
| medaka_version | String | Version of Medaka used |
| medaka_model | String | Medaka model used for consensus generation |
| samtools_version | String | Version of Samtools used |

## References

> Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191

> Medaka: Nanoporetech's tool for consensus sequence generation. https://github.com/nanoporetech/medaka