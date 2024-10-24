# Dorado Basecalling Workflow - Version 1.0

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Standalone](../../workflows_overview/workflows_type.md/#standalone) | [Any Taxa](../../workflows_overview/workflows_kingdom.md/#any-taxa) | Dorado v1.0 | Yes | Sample-level |

## Dorado Basecalling Overview

The Dorado Basecalling workflow is used to convert Oxford Nanopore `POD5` sequencing files into `FASTQ` format by utilizing a GPU-accelerated environment. This workflow is ideal for high-throughput applications where fast and accurate basecalling is essential. The workflow will upload fastq files to a user designated terra table for downstream analysis.

### Model Type Selection

Users can choose between automatic or manual model selection using a configurable use_auto_model flag:

Automatic Model Selection: Automatically picks the best model ('sup', 'hac', or 'fast') based on the input file and user-defined model accuracy paramater.

Manual Model Input: If the user disables automatic selection, a specific model path or model version must be provided.

- **Model Type (sup):** (super accuracy) The most accurate model, recommended for critical applications requiring the highest basecall accuracy. It is the slowest of the three model types.
- **Model Type (hac):** (High Accuracy) A balance between speed and accuracy, recommended for most users. Provides accurate results faster than `sup` but less accurate than `sup`.
- **Model Type (fast):** (Fast Model) The fastest model, recommended when speed is prioritized over accuracy, such as for initial analyses or non-critical applications.

### Example Manual Models:
- `dna_r10.4.1_e8.2_400bps_sup@v4.2.0`
- `dna_r10.4.1_e8.2_400bps_hac@v4.2.0`
- `dna_r10.4.1_e8.2_400bps_fast@v4.2.0`

### Workflow Structure

1. **Dorado Basecalling**: Converts `POD5` files to 'SAM' files using the specified model.
2. **Samtools Convert**: Converts the generated SAM files to BAM for efficient processing.
3. **Dorado Demultiplexing**: Demultiplexes BAM files to produce barcode-specific FASTQ files.
4. **FASTQ File Transfer**: Transfers files to Terra for downstream analysis.
5. **Terra Table Creation**: Generates a Terra table with the uploaded FASTQ files for downstream analyses.

---

## Inputs

| **Task** | **Variable** | **Type** | **Description** | **Default Value** | **Required** |
|---|---|---|---|---|---|
| Basecalling | **input_files** | Array[File] | Array of `POD5` files for basecalling | None | Yes |
| Basecalling | **use_auto_model** | Boolean | Use automatic model selection (`sup`, `hac`, or `fast` based on model accuracy)| true | No |
| Basecalling | **model_accuracy** | String | Desired model accuracy (`sup`, `hac`, `fast`) if using automatic selection | sup | No |
| Basecalling | **fastq_file_name** | String | Prefix for naming output FASTQ files | None | Yes |
| Basecalling | **dorado_model** | String | Model type (e.g., `dna_r10.4.1_e8.2_260bps_sup@v3.5.2`) if manual input | None | Yes |
| Basecalling | **kit_name** | String | Sequencing kit name used (e.g., `SQK-RPB114-24`). | None | Yes |
| Basecalling | **cpu** | Int | Number of CPUs allocated | 8 | No |
| Basecalling | **memory** | String | Amount of memory to allocate | 32GB | No |
| Basecalling | **gpuCount** | Int | Number of GPUs to use | 1 | No |
| Basecalling | **gpuType** | String | Type of GPU (e.g., `nvidia-tesla-t4`). | nvidia-tesla-t4 | No |
| Demultiplexing | **fastq_upload_path** | String | Location to upload FASTQ files on Terra (copy path from terra folder) | None | Yes |
| Demultiplexing | **fastq_file_name** | String | Prefix for naming output FASTQ files| None| Yes |
| Terra Table | **terra_project** | String | Terra project ID for final fastq file uplaod to terra table | None | Yes |
| Terra Table | **terra_workspace** | String | Terra workspace name for final fastq file upload to Terra table | None | Yes |

---

### Detailed Input Information
- **fastq_file_name**: This will serve as a prefix for the output FASTQ files. For example, if you provide `project001`, the resulting files will be named `project001_barcodeXX.fastq.gz`.
- **kit_name**: Ensure the correct kit name is provided, as it determines the barcoding and adapter trimming behavior.
- **fastq_upload_path**: This is the folder path in Terra where the final FASTQ files will be transferred for further analysis. Ensure the path matches your Terra workspace configuration.

---

## Outputs

| **Variable** | **Type** | **Description** |
|---|---|---|
| **basecalled_fastqs** | Array[File] | Array of FASTQ files generated from basecalling |
| **demuxed_fastqs** | Array[File] | FASTQ files produced from BAM demultiplexing |
| **logs** | Array[File] | Log files from the demultiplexing process |
| **terra_table_tsv** | File | TSV file for Terra table upload |

<!-- -->
><https://github.com/nanoporetech/dorado/>
