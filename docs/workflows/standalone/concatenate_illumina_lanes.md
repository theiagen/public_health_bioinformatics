# Concatenate Illumina Lanes

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Standalone](../../workflows_overview/workflows_type.md/#standalone) | [Any Taxa](../../workflows_overview/workflows_kingdom.md/#any-taxa) | PHB 2.3.0 | Yes | Sample-level |

## Concatenate_Illumina_Lanes_PHB

Some Illumina machines produce multi-lane FASTQ files for a single sample. This workflow concatenates the multiple lanes into a single FASTQ file per read type (forward or reverse).

### Inputs

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| concatenate_illumina_lanes | **read1_lane1** | File | The first lane for the forward reads | | Required |
| concatenate_illumina_lanes | **read1_lane2** | File | The second lane for the forward reads | | Required |
| concatenate_illumina_lanes | **samplename** | String | The name of the sample, used to name the output files | | Required |
| cat_lanes | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| cat_lanes | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 50 | Optional |
| cat_lanes | **docker** | String | The Docker container to use for the task |  "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.2" | Optional |
| cat_lanes | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 4 | Optional |
| concatenate_illumina_lanes | **read1_lane3** | File | The third lane for the forward reads | | Optional |
| concatenate_illumina_lanes | **read1_lane4** | File | The fourth lane for the forward reads | | Optional |
| concatenate_illumina_lanes | **read2_lane1** | File | The first lane for the reverse reads | | Optional |
| concatenate_illumina_lanes | **read2_lane2** | File | The second lane for the reverse reads | | Optional |
| concatenate_illumina_lanes | **read2_lane3** | File | The third lane for the reverse reads | | Optional |
| concatenate_illumina_lanes | **read2_lane4** | File | The fourth lane for the reverse reads | | Optional |
| version_capture | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) | | Optional |

### Workflow Tasks

This workflow concatenates the Illumina lanes for forward and (if provided) reverse reads. The output files are named as followed:

- Forward reads: `<samplename>_merged_R1.fastq.gz`
- Reverse reads: `<samplename>_merged_R2.fastq.gz`

### Outputs

| **Variable** | **Type** | **Description** |
|---|---|---|
| concatenate_illumina_lanes_analysis_date | String | Date of analysis |
| concatenate_illumina_lanes_version | String | Version of PHB used for the analysis |
| read1_concatenated | File | Concatenated forward reads |
| read2_concatenated | File | Concatenated reverse reads |
