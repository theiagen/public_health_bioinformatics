# Pangolin_Update

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Genomic Characterization](../../workflows_overview/workflows_type.md/#genomic-characterization) | [Viral](../../workflows_overview/workflows_kingdom.md/#viral), SARS-Cov-2 | PHB v2.0.0 | Yes | Sample-level |

## Pangolin_Update_PHB

The Pangolin_Update workflow re-runs Pangolin updating prior lineage calls from one docker image to meet the lineage calls specified in an alternative docker image. The most common use case for this is updating lineage calls to be up-to-date with the latest Pangolin nomenclature by using the latest available Pangolin docker image ([found here](https://www.notion.so/theiagen/Docker-Image-and-Reference-Materials-for-SARS-CoV-2-Genomic-Characterization-98328c61f5cb4f77975f512b55d09108)).

### Inputs

This workflow runs on the sample level.

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| pangolin_update | **assembly_fasta** | File | SARS-CoV-2 assembly file in FASTA format |  | Required |
| pangolin_update | **old_lineage** | String | The Pangolin lineage previously assigned to the sample |  | Required |
| pangolin_update | **old_pangolin_assignment_version** | String | Version of the Pangolin software previously used for lineage assignment. |  | Required |
| pangolin_update | **old_pangolin_docker** | String | The Pangolin docker image previously used for lineage assignment. |  | Required |
| pangolin_update | **old_pangolin_versions** | String | All pangolin software and database versions previously used for lineage assignment. |  | Required |
| pangolin_update | **samplename** | String | The name of the sample being analyzed. |  | Required |
| pangolin_update | **lineage_log** | File | TSV file detailing previous lineage assignments and software versions for this sample.  |  | Optional |
| pangolin_update | **new_pangolin_docker** | String | The Pangolin docker image used to update the Pangolin lineage assignments. |  | Optional |
| pangolin4 | **analysis_mode** | String | Pangolin inference engine for lineage designations (usher or pangolearn) | None | Optional |
| pangolin4 | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional |
| pangolin4 | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| pangolin4 | **expanded_lineage** | Boolean | True/False that determines if a lineage should be expanded without aliases (e.g., BA.1 → B.1.1.529.1) | TRUE | Optional |
| pangolin4 | **max_ambig** | Float | Maximum proportion of Ns allowed for Pangolin to attempt assignment | 0.5 | Optional |
| pangolin4 | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| pangolin4 | **min_length** | Int | Minimum query length allowed for pangolin to attempt assignment | 10000 | Optional |
| pangolin4 | **pangolin_arguments** | String | Optional arguments for pangolin e.g. "--skip-scorpio" | None | Optional |
| pangolin4 | **skip_designation_cache** | Boolean | True/False that determines if the designation cache should be used | FALSE | Optional |
| pangolin4 | **skip_scorpio** | Boolean | True/False that determines if scorpio should be skipped. | FALSE | Optional |
| pangolin_update_log | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional |
| pangolin_update_log | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| pangolin_update_log | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1 | Optional |
| pangolin_update_log | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| pangolin_update_log | **timezone** | String |  |  | Optional |
| version_capture | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional |

### Outputs

| **Variable** | **Type** | **Description** |
|---|---|---|
| **pango_lineage** | String | Pango lineage as determined by Pangolin |
| **pango_lineage_expanded** | String | Pango lineage without use of aliases; e.g., BA.1 → B.1.1.529.1 |
| **pango_lineage_log** | File | TSV file listing Pangolin lineage assignments and software versions for this sample |
| **pango_lineage_report** | File | Full Pango lineage report generated by Pangolin |
| **pangolin_assignment_version** | String | Version of the Pangolin software (e.g. PANGO or PUSHER) used for lineage assignment |
| **pangolin_conflict** | String | Number of lineage conflicts as determined by Pangolin |
| **pangolin_docker** | String | The Docker container to use for the task |
| **pangolin_notes** | String | Lineage notes as determined by Pangolin |
| **pangolin_update_analysis_date** | String | Date of analysis |
| **pangolin_update_version** | String | Version of the Public Health Bioinformatics (PHB) repository used |
| **pangolin_updates** | String | Result of Pangolin Update (lineage changed versus unchanged) with lineage assignment and date of analysis |
| **pangolin_versions** | String | All Pangolin software and database versions |