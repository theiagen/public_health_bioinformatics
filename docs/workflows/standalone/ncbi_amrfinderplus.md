# NCBI-AMRFinderPlus

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Name": "[**NCBI-AMRFinderPlus**](../workflows/standalone/ncbi_amrfinderplus.md)"}, columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level", "Dockstore"]) }}

## NCBIAMRFinderPlus_PHB

AMRFinderPlus identifies acquired antimicrobial resistance (AMR) genes, virulence genes, and stress genes.  Such AMR genes confer resistance to antibiotics, metals, biocides, heat, or acid. For some taxa (see [here](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus#--organism-option)), AMRFinderPlus will provide taxa-specific results including filtering out genes that are almost ubiquitous in the taxa (intrinsic genes) and identifying resistance-associated point mutations.  In TheiaProk, the taxon used by AMRFinderPlus is specified based on the `gambit_predicted_taxon` or a user-provided `expected_taxon`.

You can check if a gene or point mutation is in the AMRFinderPlus database [here](https://www.ncbi.nlm.nih.gov/pathogens/refgene/#), find the sequences of reference genes [here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA313047), and search the query Hidden Markov Models (HMMs) used by AMRFinderPlus to identify AMR genes and some stress and virulence proteins ([here](https://www.ncbi.nlm.nih.gov/pathogens/hmm/)). The AMRFinderPlus database is updated frequently. You can ensure you are using the most up-to-date version by specifying the docker image as a workflow input. You might like to save this docker image as a workspace data element to make this easier.

### 📋 Use Cases

- To run ONLY AMRFinderPlus software instead of running the entire TheiaProk workflow. This workflow will run much faster than the TheiaProk workflows.
- To update AMRFinderPlus results when a new version of the software and/or its database are released by the NCBI developers.

### Inputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "AMRFinderPlus"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "AMRFinderPlus"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///

## References

>Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG, Haendiges J, Haft DH, Hoffmann M, Pettengill JB, Prasad AB, Tillman GE, Tyson GH, Klimke W. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci Rep. 2021 Jun 16;11(1):12728. doi: 10.1038/s41598-021-91456-0. PMID: 34135355; PMCID: PMC8208984. <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8208984/>
<!-- -->
><https://github.com/ncbi/amr>
