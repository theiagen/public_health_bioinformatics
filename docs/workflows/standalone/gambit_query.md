# GAMBIT_Query

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line compatibility** | **Workflow type** |
|---|---|---|---|---|
| [Standalone](../../workflows_overview/workflows-type.md/#standalone) | [Bacteria](../../workflows_overview/workflows-kingdom.md/#bacteria), [Mycotics](../../workflows_overview/workflows-kingdom.md#mycotics) | PHB v2.2.0 | Yes | Sample-level |

## GAMBIT_Query_PHB

The GAMBIT_Query_PHB workflow performs taxon assignment of a genome assembly using the GAMBIT task.

### Inputs

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default attribute** | **Status** |
|---|---|---|---|---|---|
| gambit_query | assembly_fasta | File | Assembly file in FASTA format |  | Required |
| gambit_query | samplename | String | Sample name |  | Required |
| gambit | cpu | Int | Number of CPUs | 8 | Optional |
| gambit | disk_size | Int | Disk size in GB | 100 | Optional |
| gambit | memory | Int | Memory in GB | 16 | Optional |
| gambit | docker | String | Docker image to use | "us-docker.pkg.dev/general-theiagen/staphb/gambit:1.0.0" | Optional |
| gambit | gambit_db_genomes | File | Database of metadata for assembled query genomes; requires complementary signatures file. If not provided, uses default database "/gambit-db" | "gs://gambit-databases-rp/2.0.0/gambit-metadata-2.0.0-20240628.gdb" | Optional |
| gambit | gambit_db_signatures | File | Signatures file; requires complementary genomes file. If not specified, the file from the docker container will be used. | "gs://gambit-databases-rp/2.0.0/gambit-signatures-2.0.0-20240628.gs" | Optional |

### Workflow Tasks

[`GAMBIT`](https://github.com/jlumpe/gambit) determines the taxon of the genome assembly using a k-mer based approach to match the assembly sequence to the closest complete genome in a database, thereby predicting its identity. Sometimes, GAMBIT can confidently designate the organism to the species level. Other times, it is more conservative and assigns it to a higher taxonomic rank.

For additional details regarding the GAMBIT tool and a list of available GAMBIT databases for analysis, please consult the [GAMBIT](https://www.notion.so/GAMBIT-7c1376b861d0486abfbc316480046bdc?pvs=21) tool documentation.

!!! techdetails "GAMBIT Technical Details"

    |  | Links |
    | --- | --- |
    | Task | [task_gambit.wdl](https://github.com/theiagen/public_health_bacterial_genomics/blob/main/tasks/taxon_id/task_gambit.wdl) |
    | Software source code | [GAMBIT on GitHub](https://github.com/jlumpe/gambit) |
    | Software documentation | [GAMBIT ReadTheDocs](https://gambit-genomics.readthedocs.io/en/latest/) |
    | Original publication | [GAMBIT (Genomic Approximation Method for Bacterial Identification and Tracking): A methodology to rapidly leverage whole genome sequencing of bacterial isolates for clinical identification](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0277575) |

### Outputs

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

## References

> GAMBIT (Genomic Approximation Method for Bacterial Identification and Tracking): A methodology to rapidly leverage whole genome sequencing of bacterial isolates for clinical identification. Lumpe et al. PLOS ONE, 2022. DOI: [10.1371/journal.pone.0277575](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0277575)