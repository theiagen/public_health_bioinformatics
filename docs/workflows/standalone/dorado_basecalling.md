# Dorado Basecalling Workflow - Version 1.0

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Standalone](../../workflows_overview/workflows_type.md/#standalone) | [Any Taxa](../../workflows_overview/workflows_kingdom.md/#any-taxa) | Dorado v1.0 | Yes | Sample-level |

## Dorado Basecalling Overview

The Dorado Basecalling workflow is used to convert Oxford Nanopore `POD5` sequencing files into `FASTQ` format by utilizing a GPU-accelerated environment. This workflow is ideal for high-throughput applications where fast and accurate basecalling is essential.

### Inputs

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** | **Workflow** |
|---|---|---|---|---|---|---|
| basecall | **input_files** | Array[File] | Array of `POD5` files to be basecalled | None | Required | Dorado |
| basecall | **sample_names** | Array[String] | Array of sample names corresponding to the input files | None | Required | Dorado |
| basecall | **dorado_model** | String | Dorado basecalling model (e.g., `dna_r10.4.1_e8.2_260bps_sup@v3.5.2`) | None | Required | Dorado |
| basecall | **output_prefix** | String | Prefix to apply to output files | None | Required | Dorado |
| basecall | **cpu** | Int | Number of CPUs to allocate to the task | 8 | Optional | Dorado |
| basecall | **memory** | String | Amount of memory/RAM to allocate to the task | 32GB | Optional | Dorado |
| basecall | **docker** | String | The Docker container used for this task | us-docker.pkg.dev/general-theiagen/staphb/dorado:0.8.0 | Optional | Dorado |
| basecall | **gpuCount** | Int | Number of GPUs to use for basecalling | 1 | Optional | Dorado |
| basecall | **gpuType** | String | Type of GPU to use | nvidia-tesla-t4 | Optional | Dorado |

### Outputs

| **Variable** | **Type** | **Description** | **Workflow** |
|---|---|---|---|
| basecalled_fastqs | Array[File] | Array of FASTQ files generated from basecalling | Dorado |
| logs | Array[File] | Array of log files capturing the basecalling process | Dorado |
