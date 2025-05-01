# SRA_Fetch

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Data Import](../../workflows_overview/workflows_type.md/#data-import) | [Any taxa](../../workflows_overview/workflows_kingdom.md/#any-taxa) | PHB v2.2.0 | Yes | Sample-level |

## SRA_Fetch_PHB

The `SRA_Fetch` workflow downloads sequence data from NCBI's Sequence Read Archive (SRA). It requires an SRA run accession then populates the associated read files to a Terra data table.

Read files associated with the SRA run accession provided as input are copied to a Terra-accessible Google bucket. Hyperlinks to those files are shown in the "read1" and "read2" columns of the Terra data table.

### Inputs

The only required input for the SRA_Fetch workflow is an SRA run accession beginning "SRR", an ENA run accession beginning "ERR", or a DRA run accession which beginning "DRR".

Please see the [NCBI Metadata and Submission Overview](https://www.ncbi.nlm.nih.gov/sra/docs/submitmeta/) for assistance with identifying accessions. Briefly, NCBI-accessioned objects have the following naming scheme:

| STUDY | SRP# |
| --- | --- |
| SAMPLE | SRS# |
| EXPERIMENT | SRX# |
| RUN  | SRR# |

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/input_tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="SRA_Fetch", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Outputs

Read data are available either with full base quality scores (**SRA Normalized Format**) or with simplified quality scores (**SRA Lite**). The **SRA Normalized Format** includes full, per-base quality scores, whereas **base quality scores** **have been simplified in SRA Lite files.** This means that all quality scores have been artificially set to Q-30 or Q3. More information about these files can be found [here](https://www.ncbi.nlm.nih.gov/sra/docs/sra-data-formats/).

Given the lack of usefulness of SRA Lite formatted FASTQ files, we try to avoid these by preferentially searching SRA directly (SRA-Lite is more probably to be the file synced to other repositories), but sometimes downloading these files is unavoidable. To make the user aware of this, a warning column is present that is populated when an SRA-Lite file is detected.

<div class="searchable-table" markdown="1">

| **Variable** | **Type** | **Description** |
|---|---|---|
| sra_fetch_version | String | The version of the repository the SRA_Fetch workflow is in |
| sra_fetch_analysis_date | String | The date the workflow was run |
| read1 | File | File containing the forward reads |
| read2 | File | File containing the reverse reads (not available for single-end or ONT data) |
| fastq_dl_date | String | The date of the read data download |
| fastq_dl_docker | String | The docker used |
| fastq_dl_metadata | File | File containing metadata of the provided accession such as submission_accession, library_selection, instrument_platform, among others |
| fastq_dl_version | String | The version of fastq-dl used |
| fastq_dl_warning | String |  This warning field is populated if SRA-Lite files are detected. These files contain all quality encoding as Phred-30 or Phred-3. |

</div>

## References

> **fastq-dl**: Petit III, R. A., Hall, M. B., Tonkin-Hill, G., Zhu, J., & Read, T. D. fastq-dl: efficiently download FASTQ files from SRA or ENA repositories (Version 2.0.2) [Computer software]. <https://github.com/rpetit3/fastq-dl>
