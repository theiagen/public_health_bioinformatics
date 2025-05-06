# GAMBIT_Query

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Standalone](../../workflows_overview/workflows_type.md/#standalone) | [Bacteria](../../workflows_overview/workflows_kingdom.md/#bacteria), [Mycotics](../../workflows_overview/workflows_kingdom.md#mycotics) | PHB v2.2.0 | Yes | Sample-level |

## GAMBIT_Query_PHB

The GAMBIT_Query_PHB workflow performs taxon assignment of a genome assembly using the GAMBIT task.

### Inputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="Gambit_Query", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Workflow Tasks

[`GAMBIT`](https://github.com/jlumpe/gambit) determines the taxon of the genome assembly using a k-mer based approach to match the assembly sequence to the closest complete genome in a database, thereby predicting its identity. Sometimes, GAMBIT can confidently designate the organism to the species level. Other times, it is more conservative and assigns it to a higher taxonomic rank.

For additional details regarding the GAMBIT tool and a list of available GAMBIT databases for analysis, please consult the [GAMBIT](../../guides/gambit.md) tool documentation.

!!! techdetails "GAMBIT Technical Details"

    |  | Links |
    | --- | --- |
    | Task | [task_gambit.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/task_gambit.wdl) |
    | Software Source Code | [GAMBIT on GitHub](https://github.com/jlumpe/gambit) |
    | Software Documentation | [GAMBIT ReadTheDocs](https://gambit-genomics.readthedocs.io/en/latest/) |
    | Original Publication(s) | [GAMBIT (Genomic Approximation Method for Bacterial Identification and Tracking): A methodology to rapidly leverage whole genome sequencing of bacterial isolates for clinical identification](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0277575) |

### Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filter_column="Workflow", filter_values="Gambit_Query", columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///

> GAMBIT (Genomic Approximation Method for Bacterial Identification and Tracking): A methodology to rapidly leverage whole genome sequencing of bacterial isolates for clinical identification. Lumpe et al. PLOS ONE, 2022. DOI: [10.1371/journal.pone.0277575](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0277575)
