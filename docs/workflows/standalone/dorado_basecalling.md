# Dorado Basecalling Workflow - Version 1.0

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Standalone](../../workflows_overview/workflows_type.md/#standalone) | [Any Taxa](../../workflows_overview/workflows_kingdom.md/#any-taxa) | Dorado v1.0 | Yes | Sample-level |

## Dorado Basecalling Overview

The Dorado Basecalling workflow is used to convert Oxford Nanopore `POD5` sequencing files into `FASTQ` format by utilizing a GPU-accelerated environment. This workflow is ideal for high-throughput applications where fast and accurate basecalling is essential.

### Model Type Selection

- **Model Type (sup):** (super accuracy) The most accurate model, recommended for critical applications requiring the highest basecall accuracy. It is the slowest of the three model types.
- **Model Type (hac):** (High Accuracy) A balance between speed and accuracy, recommended for most users. Provides accurate results faster than `sup` but less accurate than `sup`.
- **Model Type (fast):** (Fast Model) The fastest model, recommended when speed is prioritized over accuracy, such as for initial analyses or non-critical applications.

**Example Models:**
- `dna_r10.4.1_e8.2_400bps_sup@v4.2.0`
- `dna_r10.4.1_e8.2_400bps_hac@v4.2.0`
- `dna_r10.4.1_e8.2_400bps_fast@v4.2.0`

## **Workflow Structure**

1. **Dorado Basecalling**: Converts `POD5` files to **SAM** files using the specified model.
2. **Samtools Convert**: Converts the generated SAM files to BAM for efficient processing.
3. **Dorado Demultiplexing**: Demultiplexes BAM files to produce barcode-specific FASTQ files.
4. **FASTQ File Transfer**: Transfers files to Terra for downstream analysis.
5. **Terra Table Creation**: Generates a Terra table with the uploaded FASTQ files.

---

## **Inputs**

| **Task** | **Variable** | **Type** | **Description** | **Default Value** | **Required** | **Workflow** |
|---|---|---|---|---|---|---|
| Basecalling | **input_files** | Array[File] | Array of `POD5` files for basecalling | None | Yes | Dorado |
| Basecalling | **fastq_file_name** | String | Prefix for naming output FASTQ files | None | Yes | Dorado |
| Basecalling | **dorado_model** | String | Model type (e.g., `dna_r10.4.1_e8.2_260bps_sup@v3.5.2`) | None | Yes | Dorado |
| Basecalling | **kit_name** | String | Kit name used for sequencing | None | Yes | Dorado |
| Basecalling | **cpu** | Int | Number of CPUs allocated | 8 | No | Dorado |
| Basecalling | **memory** | String | Amount of memory to allocate | 32GB | No | Dorado |
| Basecalling | **gpuCount** | Int | Number of GPUs to use | 1 | No | Dorado |
| Basecalling | **gpuType** | String | Type of GPU | nvidia-tesla-t4 | No | Dorado |
| Demultiplexing | **fastq_upload_path** | String | Location to upload FASTQ files on Terra (copy path from terra folder) | None | Yes | Demux |
| Terra Table | **terra_project** | String | Terra project ID | None | Yes | Terra |
| Terra Table | **terra_workspace** | String | Terra workspace name | None | Yes | Terra |

---

## **Outputs**

| **Variable** | **Type** | **Description** | **Workflow** |
|---|---|---|---|
| **basecalled_fastqs** | Array[File] | Array of FASTQ files generated from basecalling | Dorado |
| **demuxed_fastqs** | Array[File] | FASTQ files produced from BAM demultiplexing | Demux |
| **logs** | Array[File] | Log files from the demultiplexing process | Demux |
| **terra_table_tsv** | File | TSV file for Terra table upload | Terra |
