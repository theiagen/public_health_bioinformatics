# GAMBIT_Query

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Standalone](../../workflows_overview/workflows_type.md/#standalone) | [Bacteria](../../workflows_overview/workflows_kingdom.md/#bacteria), [Mycotics](../../workflows_overview/workflows_kingdom.md#mycotics) | PHB v2.2.0 | Yes | Sample-level |

## GAMBIT_Query_PHB

The GAMBIT_Query_PHB workflow performs taxon assignment of a genome assembly using the GAMBIT task.

### Inputs

<div class="searchable-table" markdown="1">

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| gambit_query | **assembly_fasta** | File | Assembly file in FASTA format |  | Required |
| gambit_query | **samplename** | String | Sample name |  | Required |
| gambit | **cpu** | Int | Number of CPUs to allocate to the task | 8 | Optional |
| gambit | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| gambit | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 16 | Optional |
| gambit | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/staphb/gambit:1.0.0" | Optional |
| gambit | **gambit_db_genomes** | File | Database of metadata for assembled query genomes; requires complementary signatures file. If not provided, uses default database "/gambit-db" | "gs://gambit-databases-rp/2.0.0/gambit-metadata-2.0.0-20240628.gdb" | Optional |
| gambit | **gambit_db_signatures** | File | Signatures file; requires complementary genomes file. If not specified, the file from the docker container will be used. | "gs://gambit-databases-rp/2.0.0/gambit-signatures-2.0.0-20240628.gs" | Optional |

</div>

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

<div class="searchable-table" markdown="1">

| **Variable** | **Type** | **Description** |
|---|---|---|
| gambit_closest_genomes | File | CSV file listing genomes in the GAMBIT database that are most similar to the query assembly |
| gambit_db_version | String | Version of the GAMBIT database used |
| gambit_docker | String | GAMBIT Docker used |
| gambit_predicted_taxon | String | Taxon predicted by GAMBIT |
| gambit_predicted_taxon_rank | String | Taxon rank of GAMBIT taxon prediction |
| gambit_query_wf_analysis_date | String | Date of analysis |
| gambit_query_wf_version | String | PHB repository version |
| gambit_report | File | GAMBIT report in a machine-readable format |
| gambit_version | String | Version of gambit software used |

</div>

> GAMBIT (Genomic Approximation Method for Bacterial Identification and Tracking): A methodology to rapidly leverage whole genome sequencing of bacterial isolates for clinical identification. Lumpe et al. PLOS ONE, 2022. DOI: [10.1371/journal.pone.0277575](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0277575)
