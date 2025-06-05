# Host Decontaminate

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Name": "[**Host_Decontaminate**](../workflows/standalone/host_decontaminate.md)"}, columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level"]) }}

## Host_Decontaminate_PHB

Host genetic data is frequently incidentally sequenced alongside pathogens, which can negatively affect the quality of downstream analysis. Host Decontaminate attempts to remove host reads by aligning to a reference host genome acquired on-the-fly. The reference host genome can be acquired via [NCBI Taxonomy-compatible](https://www.ncbi.nlm.nih.gov/taxonomy) taxon input or assembly accession. Host Decontaminate maps inputted reads to the host genome using `minimap2`, reports mapping statistics to this host genome, and outputs the unaligned dehosted reads. 

A common alternative approach to decontamination is to use a fast, usually k-mer-based, read classification software such as `kraken2`/`Metabuli`. These software have the advantage of screening large databases of many organisms much faster than full read alignment, though the quality of decontamination is confounded by high rates of false negative/false positive alignments compared to direct alignment. When the host is known, direct alignment to a reference host genome via `minimap2` is likely to yield higher quality mapping, which is the underlying justification for Host Decontaminate's approach.

### Inputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "Host_Decontaminate"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

??? warning "`host`"
    The `host` input needs to be an [NCBI Taxonomy-compatible](https://www.ncbi.nlm.nih.gov/taxonomy) taxon or accession. If using a taxon, the first retrieved genome corresponding to that taxon is retrieved. If using an accession, it must be coupled with the `is_accession` boolean populated as "true".

### Workflow Tasks

{{ include_md("common_text/ncbi_datasets_task.md", condition="host_decontaminate") }}

{{ include_md("common_text/parse_mapping_task.md", condition="bam_to_unaligned_fastq") }}

{{ include_md("common_text/assembly_metrics_task.md", condition="host_decontaminate", replacements={'`assembly_metrics`' : '`read_mapping_stats`'}) }}


### Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "Host_Decontaminate"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///
