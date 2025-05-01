# NCBI-AMRFinderPlus

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Standalone](../../workflows_overview/workflows_type.md/#standalone) | [Bacteria](../../workflows_overview/workflows_kingdom.md/#bacteria), [Mycotics](../../workflows_overview/workflows_kingdom.md#mycotics) | PHB v3.0.1 | Yes | Sample-level |

## NCBIAMRFinderPlus_PHB

AMRFinderPlus identifies acquired antimicrobial resistance (AMR) genes, virulence genes, and stress genes.  Such AMR genes confer resistance to antibiotics, metals, biocides, heat, or acid. For some taxa (see [here](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus#--organism-option)), AMRFinderPlus will provide taxa-specific results including filtering out genes that are almost ubiquitous in the taxa (intrinsic genes) and identifying resistance-associated point mutations.  In TheiaProk, the taxon used by AMRFinderPlus is specified based on the `gambit_predicted_taxon` or a user-provided `expected_taxon`.

You can check if a gene or point mutation is in the AMRFinderPlus database [here](https://www.ncbi.nlm.nih.gov/pathogens/refgene/#), find the sequences of reference genes [here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA313047), and search the query Hidden Markov Models (HMMs) used by AMRFinderPlus to identify AMR genes and some stress and virulence proteins ([here](https://www.ncbi.nlm.nih.gov/pathogens/hmm/)). The AMRFinderPlus database is updated frequently. You can ensure you are using the most up-to-date version by specifying the docker image as a workflow input. You might like to save this docker image as a workspace data element to make this easier.

### 📋 Use Cases

- To run ONLY AMRFinderPlus software instead of running the entire TheiaProk workflow. This workflow will run much faster than the TheiaProk workflows.
- To update AMRFinderPlus results when a new version of the software and/or its database are released by the NCBI developers.

### Inputs

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/input_tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="AMRFinderPlus", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Outputs

<div class="searchable-table" markdown="1">

| **Variable** | **Type** | **Description** |
|---|---|---|
| amrfinderplus_all_report | File | Output TSV file from AMRFinderPlus (described [here](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus#fields)) |
| amrfinderplus_amr_classes | String | AMRFinderPlus predictions for classes of drugs that genes found in the reads are known to confer resistance to |
| amrfinderplus_amr_core_genes | String | AMR genes identified by AMRFinderPlus where the scope is "core" |
| amrfinderplus_amr_plus_genes | String | AMR genes identified by AMRFinderPlus where the scope is "plus" |
| amrfinderplus_amr_report | File | TSV file detailing AMR genes only, from the amrfinderplus_all_report |
| amrfinderplus_amr_subclasses | String | More specificity about the drugs that genes identified in the reads confer resistance to |
| amrfinderplus_db_version | String | AMRFinderPlus database version used |
| amrfinderplus_stress_genes | String | Stress genes identified by AMRFinderPlus |
| amrfinderplus_stress_report | File | TSV file detailing stress genes only, from the amrfinderplus_all_report |
| amrfinderplus_version | String | AMRFinderPlus software version used |
| amrfinderplus_virulence_genes | String | Virulence genes identified by AMRFinderPlus |
| amrfinderplus_virulence_report | File | TSV file detailing virulence genes only, from the amrfinderplus_all_report |
| amrfinderplus_wf_analysis_date | String | Date of analysis |
| amrfinderplus_wf_version | String | Version of PHB used for the analysis |

</div>

## References

>Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG, Haendiges J, Haft DH, Hoffmann M, Pettengill JB, Prasad AB, Tillman GE, Tyson GH, Klimke W. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci Rep. 2021 Jun 16;11(1):12728. doi: 10.1038/s41598-021-91456-0. PMID: 34135355; PMCID: PMC8208984. <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8208984/>
<!-- -->
><https://github.com/ncbi/amr>
