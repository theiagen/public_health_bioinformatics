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

{{ input_table("docs/assets/input_tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="Assembly_Fetch", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Workflow Tasks

??? task "ReferenceSeeker Details (Optional)"

    ##### ReferenceSeeker {% raw %} {#referenceseeker} {% endraw %}

    `ReferenceSeeker` uses your draft assembly to identify closely related bacterial, viral, fungal, or plasmid genome assemblies in [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/).

    Databases that can be used with ReferenceSeeker are as follows:

      - archea:  `gs://theiagen-public-files-rp/terra/theiaprok-files/referenceseeker-archaea-refseq-205.v20210406.tar.gz`
      - bacterial (**default**): `gs://theiagen-public-files-rp/terra/theiaprok-files/referenceseeker-bacteria-refseq-205.v20210406.tar.gz`
      - fungi: `gs://theiagen-public-files-rp/terra/theiaprok-files/referenceseeker-fungi-refseq-205.v20210406.tar.gz`
      - plasmids: `gs://theiagen-public-files-rp/terra/theiaprok-files/referenceseeker-plasmids-refseq-205.v20210406.tar.gz`
      - viral: `gs://theiagen-public-files-rp/terra/theiaprok-files/referenceseeker-viral-refseq-205.v20210406.tar.gz`

    For ReferenceSeeker to identify a genome, it must meet user-specified thresholds for sequence coverage (`referenceseeker_conserved_dna_threshold`; default >= 0.69) and identity (`referenceseeker_ani_threshold`; default >= 0.95 ). 
    
    A list of closely related genomes is provided in `referenceseeker_tsv`. The reference genome that ranks highest according to ANI and conserved DNA values is considered the closest match and will be downloaded, with information about this provided in the `assembly_fetch_referenceseeker_top_hit_ncbi_accession` output.

    !!! techdetails "ReferenceSeeker Technical Details"

        |  | Links |
        | --- | --- |
        | Task | [task_referenceseeker.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/utilities/task_referenceseeker.wdl) |
        | Software Source Code | [ReferenceSeeker on GitHub](https://github.com/oschwengers/referenceseeker) |
        | Software Documentation | [ReferenceSeeker on GitHub](https://github.com/oschwengers/referenceseeker) |
        | Original Publication(s) | [ReferenceSeeker: rapid determination of appropriate reference genomes](https://joss.theoj.org/papers/10.21105/joss.01994) |

??? task "NCBI Datasets Details"

    ##### NCBI Datasets {% raw %} {#ncbi-datasets} {% endraw %}

    The [`NCBI Datasets`](https://www.ncbi.nlm.nih.gov/datasets/) task downloads specified assemblies from NCBI using either the [virus](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/virus-genome/) or [genome](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/genome/) (for all other genome types) package as appropriate.

    !!! techdetails "NCBI Datasets Technical Details"

        |  | Links |
        | --- | --- |
        | Task | [task_ncbi_datasets.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/data_import/task_ncbi_datasets.wdl) |
        | Software Source Code | [NCBI Datasets on GitHub](https://github.com/ncbi/datasets) |
        | Software Documentation | [NCBI Datasets Documentation on NCBI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/) |
        | Original Publication(s) | [Exploring and retrieving sequence and metadata for species across the tree of life with NCBI Datasets](https://doi.org/10.1038/s41597-024-03571-y) |

### Outputs

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/output_tables/assembly_fetch_out.tsv", input_table=False) }}

///

## References

> **ReferenceSeeker:** Schwengers O, Hain T, Chakraborty T, Goesmann A. ReferenceSeeker: rapid determination of appropriate reference genomes. J Open Source Softw. 2020 Feb 4;5(46):1994.
<!-- -->
> **NCBI Datasets:** O’Leary NA, Cox E, Holmes JB, et al. Exploring and retrieving sequence and metadata for species across the tree of life with NCBI Datasets. Sci Data 11, 732 (2024).
