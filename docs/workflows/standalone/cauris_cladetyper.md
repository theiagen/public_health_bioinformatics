# Cauris_CladeTyper

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Standalone](../../workflows_overview/workflows_type.md/#standalone) | [Mycotics](../../workflows_overview/workflows_kingdom.md#mycotics) | PHB v3.0.0 | Yes | Sample-level |

## Cauris_CladeTyper_PHB

The Cauris_CladeTyper_PHB Workflow is designed to assign the clade to _Candidozyma auris_ (also known as _Candida auris_) WGS assemblies based on their genomic sequence similarity to the five clade-specific reference files. Clade typing is essential for understanding the epidemiology and evolutionary dynamics of this emerging multidrug-resistant fungal pathogen.

### Inputs

<div class="searchable-table" markdown="1">

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| cauris_cladetyper | **assembly_fasta** | File | The input assembly file in FASTA format | | Required |
| cauris_cladetyper | **samplename** | String | The name of the sample being analyzed | | Required |
| cladetyper | **cpu** | Int | Number of CPUs to allocate to the task | 8 | Optional |
| cladetyper | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| cladetyper | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/biocontainers/hesslab-gambit:0.5.1--py37h8902056_0" | Optional |
| cladetyper | **kmer_size** | Int | The kmer size to use for generating the GAMBIT signatures file; see GAMBIT documentation for more details | 11 | Optional |
| cladetyper | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 16 | Optional |
| cladetyper | **ref_clade1** | File | The reference assembly for clade 1 | gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade1_GCA_002759435.2_Cand_auris_B8441_V2_genomic.fasta | Optional |
| cladetyper | **ref_clade1_annotated** | String | The path to the annotated reference for clade 1 | "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade1_GCA_002759435_Cauris_B8441_V2_genomic.gbff" | Optional |
| cladetyper | **ref_clade2** | File | The reference assembly for clade 2 | gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade2_GCA_003013715.2_ASM301371v2_genomic.fasta | Optional |
| cladetyper | **ref_clade2_annotated** | String | The path to the annotated reference for clade 2 | "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade2_GCA_003013715.2_ASM301371v2_genomic.gbff"| Optional |
| cladetyper | **ref_clade3** | File | The reference assembly for clade 3 | gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade3_reference.fasta | Optional |
| cladetyper | **ref_clade3_annotated** | String | The path to the annotated reference for clade 3 | "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade3_GCF_002775015.1_Cand_auris_B11221_V1_genomic.gbff" | Optional |
| cladetyper | **ref_clade4** | File | The reference assembly for clade 4 | gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade4_reference.fasta | Optional |
| cladetyper | **ref_clade4_annotated** | String | The path to the annotated reference for clade 4 | "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade4_GCA_003014415.1_Cand_auris_B11243_genomic.gbff" | Optional |
| cladetyper | **ref_clade5** | File | The reference assembly for clade 5 | gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade5_GCA_016809505.1_ASM1680950v1_genomic.fasta | Optional |
| cladetyper | **ref_clade5_annotated** | String | The path to the annotated reference for clade 5 | "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade5_GCA_016809505.1_ASM1680950v1_genomic.gbff" | Optional |
| version_capture | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional |

</div>

### Workflow Tasks

??? task "Cauris_Cladetyper"
    The Cauris_Cladetyper Workflow for _Candidozyma auris_ employs GAMBIT for taxonomic identification, comparing whole genome sequencing data against reference databases to accurately classify _Candidozyma auris_ isolates.

    A custom GAMBIT database is created using five clade-specific _Candidozyma auris_ reference genomes. Sequences undergo genomic signature comparison against this database, which then enables assignment to one of the five _Candidozyma auris_ clades (Clade I to Clade V) based on sequence similarity and phylogenetic relationships. This integrated approach ensures precise clade assignments, crucial for understanding the genetic diversity and epidemiology of _Candidozyma auris_.

    See more information on the reference information for the five clades below:

    | Clade | Genome Accession | Assembly Name | Strain | BioSample Accession |
    |---|---|---|---|---|
    | Clade I | GCA_002759435.2 | Cand_auris_B8441_V2 | B8441 | SAMN05379624 |
    | Clade II | GCA_003013715.2 | ASM301371v2 | B11220 | SAMN05379608 |
    | Clade III | GCA_002775015.1 | Cand_auris_B11221_V1 | B11221 | SAMN05379609 |
    | Clade IV | GCA_003014415.1 | Cand_auris_B11243 | B11243 | SAMN05379619 |
    | Clade V | GCA_016809505.1 | ASM1680950v1 | IFRC2087 | SAMN11570381 |

    !!! techdetails "Cauris_Cladetyper Technical Details"

        |  | Links |
        | --- | --- |
        | Task | [task_cauris_cladetyper.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/candida/task_cauris_cladetyper.wdl) |
        | Software Source Code | [GAMBIT on GitHub](https://github.com/jlumpe/gambit) |
        | Software Documentation | [GAMBIT Overview](https://theiagen.notion.site/GAMBIT-7c1376b861d0486abfbc316480046bdc?pvs=4) |
        | Original Publication(s) | [GAMBIT (Genomic Approximation Method for Bacterial Identification and Tracking): A methodology to rapidly leverage whole genome sequencing of bacterial isolates for clinical identification](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0277575) <br> [TheiaEuk: a species-agnostic bioinformatics workflow for fungal genomic characterization](https://doi.org/10.3389/fpubh.2023.1198213) |
    
### Outputs

<div class="searchable-table" markdown="1">

| **Variable** | **Type** | **Description** |
|---|---|---|
| cauris_cladetyper_wf_analysis_date | String | Date of analysis |
| cauris_cladetyper_wf_version | String | Version of PHB used for the analysis |
| cladetyper_annotated_reference | String | The annotated reference file for the identified clade, "None" if no clade was identified |
| cladetyper_clade | String | The clade assigned to the input assembly |
| cladetyper_docker_image | String | The Docker container used for the task |
| cladetyper_gambit_version | String | The version of GAMBIT used for the analysis |

</div>

## References

> Lumpe J, Gumbleton L, Gorzalski A, Libuit K, Varghese V, Lloyd T, et al. (2023) GAMBIT (Genomic Approximation Method for Bacterial Identification and Tracking): A methodology to rapidly leverage whole genome sequencing of bacterial isolates for clinical identification. PLoS ONE 18(2): e0277575. <https://doi.org/10.1371/journal.pone.0277575>
<!-- -->
> Ambrosio, Frank, Michelle Scribner, Sage Wright, James Otieno, Emma Doughty, Andrew Gorzalski, Danielle Siao, et al. 2023. "TheiaEuk: A Species-Agnostic Bioinformatics Workflow for Fungal Genomic Characterization." Frontiers in Public Health 11. <https://doi.org/10.3389/fpubh.2023.1198213>.
