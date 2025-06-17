# VADR_Update

## Quick Facts


| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Genomic Characterization](../../workflows_overview/workflows_type.md/#genomic-characterization) | [Viral](../../workflows_overview/workflows_kingdom.md/#viral) | PHB v2.2.0 | Yes | Sample-level |

## Vadr_Update_PHB

The VADR_Update workflow updates prior VADR assessments for each sample in line with the assessment criteria in an alternative docker image. This may be useful when samples have previously been subject to VADR alerts as updates to VADR assessment criteria may mean that the sample no longer raises concern about quality. The latest docker image SARS-CoV-2 for VADR can be found [here](https://theiagen.notion.site/Docker-Image-and-Reference-Materials-for-SARS-CoV-2-Genomic-Characterization-98328c61f5cb4f77975f512b55d09108).

Various models are available for many organisms. The following table provides an overview of the recommended container to be used and what options should be passed on to VADR.

| **Organism** | **docker** | **vadr_opts** | max_length |
| --- | --- | --- | --- |
| sars-cov-2 | "us-docker.pkg.dev/general-theiagen/staphb/vadr:1.6.3" | "--noseqnamemax --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn --out_allfasta" | 30000 |
| MPXV | "us-docker.pkg.dev/general-theiagen/staphb/vadr:1.6.3" | "--glsearch -s -r --nomisc --mkey mpxv --r_lowsimok --r_lowsimxd 100 --r_lowsimxl 2000 --alt_pass discontn,dupregin --out_allfasta --minimap2 --s_overhang 150" | 210000 |
| WNV | "us-docker.pkg.dev/general-theiagen/staphb/vadr:1.6.3" | "--mkey flavi --mdir /opt/vadr/vadr-models-flavi/ --nomisc --noprotid --out_allfasta" | 11000 |
| flu | "us-docker.pkg.dev/general-theiagen/staphb/vadr:1.6.3" | "--atgonly --xnocomp --nomisc --alt_fail extrant5,extrant3 --mkey flu" | 13500 |
| rsv_a | "us-docker.pkg.dev/general-theiagen/staphb/vadr:1.6.3" | "-r --mkey rsv --xnocomp" | 15500 |
| rsv_b | "us-docker.pkg.dev/general-theiagen/staphb/vadr:1.6.3" | "-r --mkey rsv --xnocomp" | 15500 |
| HAV | "us-docker.pkg.dev/general-theiagen/staphb/vadr:1.6.3-hav" | "-r -xnocomp -mkey hav.vadr" | 10500 |

### Inputs

Please note the default values are for SARS-CoV-2.

This workflow runs on the sample level.

<div class="searchable-table" markdown="1">

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| vadr_update | **assembly_length_unambiguous** | Int | Number of unambiguous basecalls within the consensus assembly |  | Required |
| vadr_update | **docker** | String | The Docker container to use for the task |  | Required |
| vadr_update | **genome_fasta** | File | Consensus genome assembly |  | Required |
| vadr | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| vadr | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| vadr | **max_length** | Int | Maximum length for the fasta-trim-terminal-ambigs.pl VADR script | 30000 | Optional |
| vadr | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| vadr | **min_length** | Int | Minimum length subsequence to possibly replace Ns for the fasta-trim-terminal-ambigs.pl VADR script | 50 | Optional |
| vadr | **skip_length** | Int | Minimum assembly length (unambiguous) to run vadr | 10000 | Optional |
| vadr | **vadr_opts** | String | Options for the v-annotate.pl VADR script | ''--glsearch -s -r --nomisc --mkey sarscov2 --alt_fail lowscore,fstukcnf,insertnn,deletinn --mdir /opt/vadr/vadr-models/'' | Optional |
| version_capture | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional |

</div>

### Outputs

| **Variable** | **Type** | **Description** |
|---|---|---|
| vadr_alerts_list | File | File containing all of the fatal alerts as determined by VADR |
| vadr_docker | String | Docker image used to run VADR |
| vadr_fastas_zip_archive | File | Archive file (in zip format) of all VADR outputs |
| vadr_num_alerts | String | Number of fatal alerts as determined by VADR |
| vadr_update_analysis_date | String | Date of analysis |
| vadr_update_version | String | Version of the Public Health Bioinformatics (PHB) repository used |
