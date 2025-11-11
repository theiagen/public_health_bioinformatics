# kSNP4

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Name": "[**kSNP4**](../workflows/phylogenetic_construction/ksnp4.md)"}, columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level", "Dockstore"]) }}

## kSNP4_PHB

The kSNP4 workflow is for phylogenetic analysis of bacterial genomes using single nucleotide polymorphisms (SNPs) and is significantly faster and more memory efficient than its predecessor, kSNP3. There are no significant algorithmic changes between the two versions, and most modifications are transparent to the user. The kSNP4 workflow identifies SNPs amongst a set of genome assemblies, then calculates a number of phylogenetic trees based on those SNPs:

- **Pan-genome phylogenetic trees:** The term "pan-genome" is used here to describe the collective genetic content amongst the set of genomes, including regions outside of genes and other coding sequences.  Outputs based on the pan-genome are labeled with `_pan`.
- **Core-genome phylogenetic trees:** The kSNP4 workflow will also generate phylogenetic trees based on the core genome (genetic content that is present in all members of the set of genomes). Outputs based on the core-genome are labeled with `_core`.

This workflow also features an optional module, `summarize_data` that creates a presence/absence matrix for the analyzed samples from a list of indicated columns (such as AMR genes, plasmid types etc.). If the `phandango_coloring` variable is set to `true`, this will be formatted for visualization in [Phandango](https://jameshadfield.github.io/phandango/#/), else it can be viewed in Excel.

While kSNP4 introduces enhancements, much of the foundational information from kSNP3 remains relevant. You can learn more about the kSNP3 workflow, including how to visualize the outputs with MicrobeTrace in the following video, which is still applicable to kSNP4: **📺 [Using kSNP3 in Terra and Visualizing Bacterial Genomic Networks in MicrobeTrace](https://www.youtube.com/watch?v=iRpNDun46R8)**

### Inputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "kSNP4"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Workflow Actions

{{ include_md("common_text/ksnp_task.md", condition="ksnp4") }}
{{ include_md("common_text/snp_dists_task.md", condition="ksnp") }}
{{ include_md("common_text/data_summary_task.md") }}

### Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "kSNP4"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///

## References

>Barry G Hall, Jeremiah Nisbet, Building Phylogenetic Trees From Genome Sequences With kSNP4, Molecular Biology and Evolution, Volume 40, Issue 11, November 2023, msad235, <https://doi.org/10.1093/molbev/msad235>
<!-- -->
<https://github.com/tseemann/snp-dists>
