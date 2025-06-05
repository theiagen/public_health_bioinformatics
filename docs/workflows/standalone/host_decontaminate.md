# Host Decontaminate

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Name": "[**Host_Decontaminate**](../workflows/standalone/host_decontaminate.md)"}, columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level"]) }}

## Host_Decontaminate_PHB

Host reads are frequently sequenced alongside pathogens and need to be removed to generate high quality *de novo* assemblies. Host Decontaminate attempts to remove host reads by aligning to a reference host genome acquired on-the-fly via an NCBI Taxonomy-compatible taxon input or assembly accession. Host Decontaminate maps these reads using `minimap2` reports mapping statistics to this host genome, and outputs the unaligned dehosted reads. 

One approach to decontamination is to use a fast, usually k-mer-based, read classification software such as `kraken2`/`Metabuli`. These software have the advantage of being able to screen large databases of many organisms much faster than full read alignment, though they can be confounded by relatively high rates of false negative/false positive alignments. When the host is known, direct alignment to a reference host genome via `minimap2` is likely to yield higher quality mapping, which is the underlying justification for Host Decontaminate's approach.

### Inputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "Host_Decontaminate"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Workflow Tasks

<div class="grid cards" markdown>

    -   {{ include_md("common_text/ncbi_datasets_task.md", condition="host_decontaminate", indent=12) }}

</div>

<div class="grid cards" markdown>

    -   {{ include_md("common_text/parse_mapping_task.md", condition="bam_to_unaligned_fastq", indent=12) }}

</div>


<div class="grid cards" markdown>

    -   {{ include_md("common_text/assembly_metrics_task.md", condition="host_decontaminate", indent=12, replacements={'`assembly_metrics`' : '`read_mapping_stats`'}) }}

</div>

### Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "Host_Decontaminate"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///
