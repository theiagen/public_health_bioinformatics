# Cauris_CladeTyper

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Name": "[**Cauris_CladeTyper**](../workflows/standalone/cauris_cladetyper.md)"}, columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level", "Dockstore"]) }}

## Cauris_CladeTyper_PHB

The Cauris_CladeTyper_PHB Workflow is designed to assign the clade to _Candidozyma auris_ (also known as _Candida auris_) WGS assemblies based on their genomic sequence similarity to the six clade-specific reference files. Clade typing is essential for understanding the epidemiology and evolutionary dynamics of this emerging multidrug-resistant fungal pathogen.

### Inputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "Cauris_CladeTyper"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Workflow Tasks

{{ include_md("common_text/cauris_cladetyper.md", condition="cauris_cladetyper") }}

### Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "Cauris_CladeTyper"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///

## References

> Lumpe J, Gumbleton L, Gorzalski A, Libuit K, Varghese V, Lloyd T, et al. (2023) GAMBIT (Genomic Approximation Method for Bacterial Identification and Tracking): A methodology to rapidly leverage whole genome sequencing of bacterial isolates for clinical identification. PLoS ONE 18(2): e0277575. <https://doi.org/10.1371/journal.pone.0277575>
<!-- -->
> Ambrosio, Frank, Michelle Scribner, Sage Wright, James Otieno, Emma Doughty, Andrew Gorzalski, Danielle Siao, et al. 2023. "TheiaEuk: A Species-Agnostic Bioinformatics Workflow for Fungal Genomic Characterization." Frontiers in Public Health 11. <https://doi.org/10.3389/fpubh.2023.1198213>.
