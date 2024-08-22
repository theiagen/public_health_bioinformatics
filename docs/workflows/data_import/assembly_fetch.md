---
hide:
 - toc
---

# Assembly Fetch

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line compatibliity** | **Workflow type** |
|---|---|---|---|---|
| [Data Import](../../../workflows_overview/workflows-type/#data-import) | [Any taxa](../../../workflows_overview/workflows-kingdom/#any-taxa) | PHB v1.3.0 | Yes | Sample-level |

## Assembly_Fetch_PHB

The `Assembly_Fetch` workflow downloads assemblies from NCBI. This is particularly useful when you need to align reads against a reference genome, for example during a reference-based phylogenetics workflow. This workflow can be run in two ways:

1. You can provide an accession for the specific assembly that you want to download, and `Assembly_Fetch` will run only the NCBI genome download task to download this assembly,
2. You can provide an assembly, and `Assembly_Fetch` will first use the `ReferenceSeeker` task to first find the closest reference genome in RefSeq to your query assembly and then will run the NCBI genome download task to download that reference assembly.

!!! tip

    NOTE: If using Assembly_Fetch workflow version 1.3.0 or higher, the call-caching feature of Terra has been DISABLED to ensure that the workflow is run from the beginning and data is downloaded fresh. Call-caching will not be enabled, even if the user checks the box ✅ in the Terra workflow interface.

### Inputs

Assembly_Fetch requires the input samplename, and either the accession for a reference genome to download (ncbi_accession) or an assembly that can be used to query RefSeq for the closest reference genome to download (assembly_fasta).

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default attribute** | **Terra Status** |
|---|---|---|---|---|---|
| reference_fetch | **samplename** | String | Your sample’s name |  | Required |
| reference_fetch | **assembly_fasta** | File | Assembly FASTA file of your sample |  | Optional |
| reference_fetch | **ncbi_accession** | String | NCBI accession passed to the NCBI datasets task to be downloaded. Example: GCF_000006945.2 (Salmonella enterica subsp. enterica, serovar Typhimurium str. LT2 reference genome) |  | Optional |
| ncbi_datasets_download_genome_accession | **cpu** | Int | number of cpus used to run ncbi datasets | 1 | Optional |
| ncbi_datasets_download_genome_accession | **disk_size** | Int | Amount of storage in GigaBytes (GB) requested for the VM to run the NCBI datasets task | 50 | Optional |
| ncbi_datasets_download_genome_accession | **docker** | String | Docker image used to run the NCBI datasets task | "us-docker.pkg.dev/general-theiagen/staphb/ncbi-datasets:14.13.2” | Optional |
| ncbi_datasets_download_genome_accession | **include_gbff** | Boolean | set to true if you would like the GenBank Flat File (GBFF) file included in the output. It contains nucleotide sequence, metadata, and annotations. | FALSE | Optional |
| ncbi_datasets_download_genome_accession | **include_gff3** | Boolean | set to true if you would like the Genomic Feature File v3 (GFF3) file included in the output. It contains nucleotide sequence, metadata, and annotations | FALSE | Optional |
| ncbi_datasets_download_genome_accession | **memory** | Int | Amount of RAM/memory requested for running the NCBI datasets task | 4 | Optional |
| referenceseeker | **cpu** | Int | number of cpus used to run referenceseeker | 4 | Optional |
| referenceseeker | **disk_size** | Int | Amount of storage in GigaBytes (GB) requested for the VM to run the referenceseeker task | 200 | Optional |
| referenceseeker | **docker** | String | Docker image used to run the referenceseeker task | "us-docker.pkg.dev/general-theiagen/biocontainers/referenceseeker:1.8.0--pyhdfd78af_0” | Optional |
| referenceseeker | **memory** | Int | Amount of RAM/memory requested to run the referenceseeker task | 16 | Optional |
| referenceseeker | **referenceseeker_ani_threshold** | Float | ANI threshold used to exclude ref genomes when ANI value less than this value. | 0.95 | Optional |
| referenceseeker | **referenceseeker_conserved_dna_threshold** | Float | Conserved DNA threshold used to exclude ref genomes when conserved DNA value is less than this value. | 0.69 | Optional |
| referenceseeker | **referenceseeker_db** | File | Database used by the referenceseeker tool that contains bacterial genomes from RefSeq release 205. Downloaded from referenceseeker GitHub repo. | "gs://theiagen-public-files-rp/terra/theiaprok-files/referenceseeker-bacteria-refseq-205.v20210406.tar.gz” | Optional |
| version_capture | **docker** | String | The Docker image used to run the version_capture task | "us-docker.pkg.dev/general-theiagen/ubuntu/ubuntu:jammy-20230816" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional |

# Tasks

## ReferenceSeeker (optional)

`ReferenceSeeker` uses your draft assembly to identify closely related bacterial, viral, fungal, or plasmid genome assemblies in [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/). 

Databases for use with ReferenceSeeker are as follows, and can be used by pasting the gs uri in double quotation marks `" "` into the `referenceseeker_db` optional input:

- archea:  `gs://theiagen-public-files-rp/terra/theiaprok-files/referenceseeker-archaea-refseq-205.v20210406.tar.gz`
- bacterial (**default**): `gs://theiagen-public-files-rp/terra/theiaprok-files/referenceseeker-bacteria-refseq-205.v20210406.tar.gz`
- fungi: `gs://theiagen-public-files-rp/terra/theiaprok-files/referenceseeker-fungi-refseq-205.v20210406.tar.gz`
- plasmids: `gs://theiagen-public-files-rp/terra/theiaprok-files/referenceseeker-plasmids-refseq-205.v20210406.tar.gz`
- viral: `gs://theiagen-public-files-rp/terra/theiaprok-files/referenceseeker-viral-refseq-205.v20210406.tar.gz`

For ReferenceSeeker to identify a genome, it must meet user-specified thresholds for sequence coverage (`referenceseeker_conserved_dna_threshold`) and identity (`referenceseeker_ani_threshold`). The default values for these are set according to community standards (conserved DNA >= 69 % and ANI >= 95 %). A list of closely related genomes is provided in `referenceseeker_tsv`. The reference genome that ranks highest according to ANI and conserved DNA values is considered the closest match and will be downloaded, with information about this provided in the `assembly_fetch_referenceseeker_top_hit_ncbi_accession` output.

- **Inputs**
    
    [Untitled Database](https://www.notion.so/7106749fff9d49e3892b06481b7de74e?pvs=21)
    
- **Outputs**
    
    [Untitled Database](https://www.notion.so/3ac0873285ec43f78a070c84d7687ba8?pvs=21)
    
- **Technical details**
    
    
    |  | Links |
    | --- | --- |
    | Task | [task_referenceseeker.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/task_referenceseeker.wdl) |
    | Software version | 1.8.0 ([quay.io/biocontainers/referenceseeker:1.8.0--pyhdfd78af_0](http://quay.io/biocontainers/referenceseeker:1.8.0--pyhdfd78af_0)) |
    | Software source code | ‣ |
    | Software documentation | ‣ |
    | Original publication | [ReferenceSeeker: rapid determination of appropriate reference genomes](https://joss.theoj.org/papers/10.21105/joss.01994) |

## NCBI datasets

The [`NCBI datasets`](https://www.ncbi.nlm.nih.gov/datasets/) task downloads specified assemblies from NCBI using either the [virus](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/virus-genome/) or [genome](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/genome/) (for all other genome types) package as appropriate. 

- **Inputs**
    
    [Untitled Database](https://www.notion.so/d519968fad6e4f3984fc99e8e44f7d95?pvs=21)
    
- **Outputs**
    
    [Untitled Database](https://www.notion.so/dcf8ce4ec4c44cd3a1f139ad01f29fc3?pvs=21)
    
- **Technical details**
    
    
    |  | Links |
    | --- | --- |
    | Task | [task_ncbi_datasets.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/task_ncbi_datasets.wdl) |
    | Software version | 14.13.2 (staphb/ncbi-datasets:14.13.2) |
    | Software source code | ‣ |
    | Software documentation | ‣ |
    | Original publication | Not known to be published |


# Outputs

[Assembly_Fetch Outputs](https://www.notion.so/60147f89a65c4b928b51e6a12904cfef?pvs=21)

# References

**ReferenceSeeker:** Schwengers O, Hain T, Chakraborty T, Goesmann A. ReferenceSeeker: rapid determination of appropriate reference genomes. J Open Source Softw. 2020 Feb 4;5(46):1994.

**NCBI datasets: datasets:** NCBI Datasets is an experimental resource for finding and building datasets [Internet]. Github; [cited 2023 Apr 19]. Available from: https://github.com/ncbi/datasets



<!-- definitions for workflow type column -->
*[Sample-level]: This workflow is run once for each sample
*[Set-level]: This workflow is run once on a group of samples

<!-- definitions for taxa column -->
*[Any taxa]: This workflow is organism-agnostic and can be run with any taxa
*[Viral]: This workflow is compatible with any viral pathogen
*[Bacteria]: This workflow is compatible with any bacterial pathogen
*[Mycotics]: This workflow is compatible with mycotic pathogens

<!-- definition for command-line compatible column -->
*[Command-line compatible]: Command-line compatibility is determined if the workflow can be run on a local command-line environment, providing all dependencies are installed, with either `miniwdl` or `cromwell`.
*[Some optional features incompatible]: Some optional features of this workflow are incompatible with command-line use and require modification
*[Yes]: This workflow is compatible with command-line use
*[No]: This workflow is not compatible with command-line use even with modifications
