# Fetch SRR Accession Workflow

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Public Data Sharing](../../workflows_overview/workflows_type.md/#public-data-sharing) | [Any Taxa](../../workflows_overview/workflows_kingdom.md/#any-taxa) | PHB v2.3.0 | Yes | Sample-level |

## Fetch SRR Accession

This workflow retrieves the Sequence Read Archive (SRA) accession (SRR) associated with a given sample accession. The primary inputs are BioSample IDs (e.g., SAMN00000000) or SRA Experiment IDs (e.g., SRX000000), which link to sequencing data in the SRA repository.

The workflow uses the fastq-dl tool to fetch metadata from SRA and specifically parses this metadata to extract the associated SRR accession and outputs the SRR accession.

### Inputs

| **Terra Task Name** | **Variable** | **Type** | **Description**| **Default Value** | **Terra Status** |
| --- | --- | --- | --- | --- | --- |
| fetch_srr_metadata | **sample_accession** | String |  SRA-compatible accession, such as a **BioSample ID** (e.g., "SAMN00000000") or **SRA Experiment ID** (e.g., "SRX000000"), used to retrieve SRR metadata. | | Required |
| fetch_srr_metadata | **cpu** | Int | Number of CPUs allocated for the task. | 2 | Optional |
| fetch_srr_metadata | **disk_size** | Int | Disk space in GB allocated for the task. | 10 | Optional |
| fetch_srr_metadata | **docker**| String | Docker image for metadata retrieval. | `us-docker.pkg.dev/general-theiagen/biocontainers/fastq-dl:2.0.4--pyhdfd78af_0` | Optional |
| fetch_srr_metadata | **memory** | Int | Memory in GB allocated for the task. | 8 | Optional |
| version_capture | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) | | Optional |

### Workflow Tasks

This workflow has a single task that performs metadata retrieval for the specified sample accession.

??? task "`fastq-dl`: Fetches SRR metadata for sample accession"
    When provided a BioSample accession or SRA experiment ID, 'fastq-dl' collects metadata and returns the appropriate SRR accession.

    !!! techdetails "fastq-dl Technical Details"
        |  | Links | 
        | --- | --- | 
        | Task | [Task on GitHub](https://github.com/theiagen-org/phb-workflows/blob/main/tasks/utilities/data_handling/task_fetch_srr_metadata.wdl) |
        | Software Source Code | [fastq-dl Source](https://github.com/rvalieris/fastq-dl) |
        | Software Documentation | [fastq-dl Documentation](https://github.com/rvalieris/fastq-dl#documentation) |
        | Original Publication | [fastq-dl: A fast and reliable tool for downloading SRA metadata](https://doi.org/10.1186/s12859-021-04346-3) |

### Outputs

| **Variable** | **Type** | **Description**|
|---|---|---|
| srr_accession| String | The SRR accession's associated with the input sample accession.|
| fetch_srr_accession_version | String | The version of the fetch_srr_accession workflow. |
| fetch_srr_accession_analysis_date | String | The date the fetch_srr_accession analysis was run. |

## References

> Valieris, R. et al., "fastq-dl: A fast and reliable tool for downloading SRA metadata." Bioinformatics, 2021.
