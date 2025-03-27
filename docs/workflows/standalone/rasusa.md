# RASUSA

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Standalone](../../workflows_overview/workflows_type.md/#standalone) | [Any Taxa](../../workflows_overview/workflows_kingdom.md/#any-taxa) | PHB v3.0.0 | Yes | Sample-level |

## RASUSA_PHB

RASUSA functions to randomly downsample the number of raw reads to a user-defined threshold.

### 📋 Use Cases

- to reduce computing resources when samples end up with drastically more data than needed to perform analyses
- to perform limit of detection (LOD) studies to identify appropriate minimum coverage thresholds required to perform downstream analyses

### 🔧 Desired size may be specified by inputting any one of the following

- coverage (e.g. 20X)
- number of bases (e.g. "5m" for 5 megabases)
- number of reads (e.g. 100000 total reads)
- fraction of reads (e.g. 0.5 samples half the reads)

!!! info "Call-caching disabled"
    If using RASUSA_PHB workflow version v2.0.0 or higher, **the call-caching feature of Terra has been DISABLED to ensure that the workflow is run from the beginning and data is downloaded fresh.** Call-caching will not be enabled, even if the user checks the box ✅ in the Terra workflow interface.

### Inputs

<div class="searchable-table" markdown="1">

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Attribute** | **Terra Status** |
|---|---|---|---|---|---|
| rasusa_workflow | **coverage** | Float | Use to specify the desired coverage of reads after downsampling; actual coverage of subsampled reads will not be exact and may be slightly higher; always check the estimated clean coverage after performing downstream workflows to verify coverage values, when necessary | | Required |
| rasusa_workflow | **genome_length** | String | Input the approximate genome size expected in quotations; this is used to determine the number of bases required to achieve the desired coverage; acceptable metric suffixes include: `b`, `k`, `m`, `g`, and `t` for base, kilo, mega, giga, and tera, respectively | | Required |
| rasusa_workflow | **read1** | File | FASTQ file containing read1 sequences | | Required |
| rasusa_workflow | **read2** | File | FASTQ file containing read2 sequences | | Required |
| rasusa_workflow | **samplename** | String | Name of the sample to be analyzed | | Required |
| rasusa_task | **bases** | String | Explicitly define the number of bases required in the downsampled reads in quotations; when used, genome size and coverage are ignored; acceptable metric suffixes include: `b`, `k`, `m`, `g`, and `t` for base, kilo, mega, giga, and tera, respectively | | Optional |
| rasusa_task | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional |
| rasusa_task | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| rasusa_task | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/staphb/rasusa:2.1.0" | Optional |
| rasusa_task | **frac** | Float | Explicitly define the fraction of reads to keep in the subsample; when used, genome size and coverage are ignored; acceptable inputs include whole numbers and decimals, e.g. 50.0 will leave 50% of the reads in the subsample | | Optional |
| rasusa_task | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| rasusa_task | **num** | Int | Optional: explicitly define the number of reads in the subsample; when used, genome size and coverage are ignored; acceptable metric suffixes include: `b`, `k`, `m`, `g`, and `t` for base, kilo, mega, giga, and tera, respectively | | Optional |
| rasusa_task | **seed** | Int | Use to assign a name to the "random seed" that is used by the subsampler; i.e. this allows the exact same subsample to be produced from the same input file/s in subsequent runs when providing the seed identifier; do not input values for random downsampling | | Optional |
| version_capture | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional |

</div>

### Outputs

<div class="searchable-table" markdown="1">

| **Variable** | **Type** | **Description** |
|---|---|---|
| rasusa_version | String | Version of RASUSA used for the analysis |
| rasusa_wf_analysis_date | String | Date of analysis |
| rasusa_wf_version | String | Version of PHB used for the analysis |
| read1_subsampled | File | New read1 FASTQ files downsampled to desired coverage |
| read2_subsampled | File | New read2 FASTQ files downsampled to desired coverage |

</div>

!!! tip "Don't Forget!"
    Remember to use the subsampled reads in downstream analyses with `this.read1_subsampled` and `this.read2_subsampled` inputs.

!!! info "Verify"
    Confirm reads were successfully subsampled before downstream analyses by comparing read file size/s to the original read file size/s

    _View file sizes by clicking on the read file listed in the Terra data table and looking at the file size_

## References

> Hall, M. B., (2022). Rasusa: Randomly subsample sequencing reads to a specified coverage. Journal of Open Source Software, 7(69), 3941, <https://doi.org/10.21105/joss.03941>
