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

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/input_tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="RASUSA", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Outputs

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/tables/all_outputs.tsv", input_table=False, filter_column="Workflow", filter_values="RASUSA", columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///

!!! tip "Don't Forget!"
    Remember to use the subsampled reads in downstream analyses with `this.read1_subsampled` and `this.read2_subsampled` inputs.

!!! info "Verify"
    Confirm reads were successfully subsampled before downstream analyses by comparing read file size/s to the original read file size/s

    _View file sizes by clicking on the read file listed in the Terra data table and looking at the file size_

## References

> Hall, M. B., (2022). Rasusa: Randomly subsample sequencing reads to a specified coverage. Journal of Open Source Software, 7(69), 3941, <https://doi.org/10.21105/joss.03941>
