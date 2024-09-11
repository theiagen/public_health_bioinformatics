# Rename_FASTQ

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Standalone](../../workflows_overview/workflows_type.md/#standalone) | [Any Taxa](../../workflows_overview/workflows_kingdom.md/#any-taxa) | PHB v2.1.0 | Yes | Sample-level |

## Rename_FASTQ_PHB

This sample-level workflow receives a read file or a pair of read files (FASTQ), compressed or uncompressed, and returns a new, renamed and compressed FASTQ file.

### Inputs

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| rename_fastq_files | **new_filename** | String | New name for the FASTQ file(s) | | Required |
| rename_fastq_files | **read1** | File | FASTQ file containing read1 sequences | | Required |
| rename_fastq_files | **read2** | File | FASTQ file containing read2 sequences | | Optional |
| rename_PE_files or rename_SE_files | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| rename_PE_files or rename_SE_files | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| rename_PE_files or rename_SE_files | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/ubuntu/ubuntu:jammy-20230816" | Optional |
| rename_PE_files or rename_SE_files | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional |
| version_capture | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional |

### Outputs

If a reverse read (`read2`) is provided, the files get renamed to the provided `new_filename` input with the notation `<new_filename>_R1.fastq.gz` and `<new_filename>_R2.fastq.gz`. If only `read1` is provided, the file is renamed to `<new_filename>.fastq.gz`. 

| **Variable** | **Type** | **Description** |
|---|---|---|
| read1_renamed | File | New read1 FASTQ file renamed to desired filename |
| read2_renamed | File | New read2 FASTQ file renamed to desired filename |
| rename_fastq_files_analysis_date | String | Date of analysis |
| rename_fastq_files_version | String | Version of PHB used for the analysis |
