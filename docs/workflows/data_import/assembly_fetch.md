# Assembly Fetch

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Data Import](../../workflows_overview/workflows_type.md/#data-import) | [Any taxa](../../workflows_overview/workflows_kingdom.md/#any-taxa) | PHB v1.3.0 | Yes | Sample-level |

## Assembly_Fetch_PHB

The `Assembly_Fetch` workflow downloads assemblies from NCBI. This is particularly useful when you need to align reads against a reference genome, for example during a reference-based phylogenetics workflow. This workflow can be run in two ways:

1. You can provide an accession for the specific assembly that you want to download, and `Assembly_Fetch` will run only the NCBI genome download task to download this assembly,
2. You can provide an assembly, and `Assembly_Fetch` will first use the `ReferenceSeeker` task to first find the closest reference genome in RefSeq to your query assembly and then will run the NCBI genome download task to download that reference assembly.

!!! info "Call-Caching Disabled"

    If using Assembly_Fetch workflow version 1.3.0 or higher, the call-caching feature of Terra has been DISABLED to ensure that the workflow is run from the beginning and data is downloaded fresh. Call-caching will not be enabled, even if the user checks the box ✅ in the Terra workflow interface.

### Inputs

Assembly_Fetch requires the input samplename, and either the accession for a reference genome to download (ncbi_accession) or an assembly that can be used to query RefSeq for the closest reference genome to download (assembly_fasta).

This workflow runs on the sample level.

!!! warning "Note on Downloading Viral Assemblies"

    If downloading viral assemblies, set `use_ncbi_virus` to true.

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="Assembly_Fetch", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Workflow Tasks

{{ include_md("common_text/referenceseeker_task.md") }}
{{ include_md("common_text/ncbi_datasets_task.md") }}

### Outputs

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/tables/all_outputs.tsv", input_table=False, filter_column="Workflow", filter_values="Assembly_Fetch", columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///

## References

> **ReferenceSeeker:** Schwengers O, Hain T, Chakraborty T, Goesmann A. ReferenceSeeker: rapid determination of appropriate reference genomes. J Open Source Softw. 2020 Feb 4;5(46):1994.
<!-- -->
> **NCBI Datasets:** O’Leary NA, Cox E, Holmes JB, et al. Exploring and retrieving sequence and metadata for species across the tree of life with NCBI Datasets. Sci Data 11, 732 (2024).
