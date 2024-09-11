# SRA_Fetch

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Data Import](../../workflows_overview/workflows_type.md/#data-import) | [Any taxa](../../workflows_overview/workflows_kingdom.md/#any-taxa) | PHB v2.2.0 | Yes | Sample-level |

## SRA_Fetch_PHB

The `SRA_Fetch` workflow downloads sequence data from NCBI's Sequence Read Archive (SRA). It requires an SRA run accession then populates the associated read files to a Terra data table.

Read files associated with the SRA run accession provided as input are copied to a Terra-accessible Google bucket. Hyperlinks to those files are shown in the "read1" and "read2" columns of the Terra data table.

### Inputs

This workflow runs on the sample level.

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| fetch_sra_to_fastq | **sra_accession** | String | SRA, ENA, or DRA accession number | | Required |
| fetch_sra_to_fastq | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| fetch_sra_to_fastq | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| fetch_sra_to_fastq | **docker_image** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/biocontainers/fastq-dl:2.0.4--pyhdfd78af_0" | Optional |
| fetch_sra_to_fastq | **fastq_dl_options** | String | Additional parameters to pass to fastq_dl from [here](https://github.com/rpetit3/fastq-dl?tab=readme-ov-file#usage) | "--provider sra" | Optional |
| fetch_sra_to_fastq | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |

The only required input for the SRA_Fetch workflow is an SRA run accession beginning "SRR", an ENA run accession beginning "ERR", or a DRA run accession which beginning "DRR".

Please see the [NCBI Metadata and Submission Overview](https://www.ncbi.nlm.nih.gov/sra/docs/submitmeta/) for assistance with identifying accessions. Briefly, NCBI-accessioned objects have the following naming scheme:

| STUDY | SRP# |
| --- | --- |
| SAMPLE | SRS# |
| EXPERIMENT | SRX# |
| RUN  | SRR# |

### Outputs

Read data are available either with full base quality scores (**SRA Normalized Format**) or with simplified quality scores (**SRA Lite**). The **SRA Normalized Format** includes full, per-base quality scores, whereas **base quality scores** **have been simplified in SRA Lite files.** This means that all quality scores have been artificially set to Q-30 or Q3. More information about these files can be found [here](https://www.ncbi.nlm.nih.gov/sra/docs/sra-data-formats/).

Given the lack of usefulness of SRA Lite formatted FASTQ files, we try to avoid these by selecting as provided SRA directly (SRA-Lite is more probably to be the file synced to other repositories), but some times downloading these files is unavoidable. To make the user aware of this, a warning column is present that is populated when an SRA-Lite file is detected.

| **Variable** | **Type** | **Description** | **Production Status** |
|---|---|---|---|
| read1 | File | File containing the forward reads | Always produced |
| read2 | File | File containing the reverse reads (not availablae for single-end or ONT data) | Produced only for paired-end data |
| fastq_dl_date | String | The date of download | Always produced |
| fastq_dl_docker | String | The docker used | Always produced |
| fastq_dl_metadata | File | File containing metadata of the provided accession such as submission_accession, library_selection, instrument_platform, among others | Always produced |
| fastq_dl_version | String | Fastq_dl version used | Always produced |
| fastq_dl_warning | String |  This warning field is populated if SRA-Lite files are detected. These files contain all quality encoding as Phred-30 or Phred-3. | Depends on internal workflow logic |

## References

> This workflow relies on [fastq-dl](https://github.com/rpetit3/fastq-dl), a very handy bioinformatics tool by Robert A. Petit III
