
??? task "Cladetyping: clade determination"
<!-- if: cauris_cladetyper -->
    The Cauris_Cladetyper Workflow for _Candidozyma auris_ employs GAMBIT for taxonomic identification, comparing whole genome sequencing data against reference databases to accurately classify _Candidozyma auris_ isolates.
<!-- endif -->

    A custom GAMBIT database is created using six clade-specific _Candidozyma auris_ reference genomes. Sequences undergo genomic signature comparison against this database, which then enables assignment to one of the six _Candidozyma auris_ clades (Clade I to Clade VI) based on sequence similarity and phylogenetic relationships. This integrated approach ensures precise clade assignments, crucial for understanding the genetic diversity and epidemiology of _Candidozyma auris_.

    See more information on the reference information for the six clades below:

    | Clade | Genome Accession | Assembly Name | Strain | BioSample Accession |
    |---|---|---|---|---|
    | Clade I | GCA_002759435.3 | Cand_auris_B8441_V3 | B8441 | SAMN05379624 |
    | Clade II | GCA_003013715.2 | ASM301371v2 | B11220 | SAMN05379608 |
    | Clade III | GCA_002775015.1 | Cand_auris_B11221_V1 | B11221 | SAMN05379609 |
    | Clade IV | GCA_003014415.1 | Cand_auris_B11243 | B11243 | SAMN05379619 |
    | Clade V | GCA_016809505.1 | ASM1680950v1 | IFRC2087 | SAMN11570381 |
    | Clade VI | GCA_032714025.1 | ASM3271402v1 | F1580 | SAMN36753179 |

    !!! warning "Clade VI annotation"
        
        Clade VI does not have an available reference genome annotation at the time of adding the reference genome into Cladetyping. While Clade VI assignment is functional, downstream variant calling is not currently possible without an annotation. Users may provide a close relative annotation, such as Clade IV, though it is unknown if Clade VI variants can reliably be called with respect to such a reference. 

    !!! techdetails "Cauris_Cladetyper Technical Details"

        |  | Links |
        | --- | --- |
        | Task | [task_cauris_cladetyper.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/candidozyma/task_cauris_cladetyper.wdl) |
        | Software Source Code | [GAMBIT on GitHub](https://github.com/jlumpe/gambit) |
        | Software Documentation | [GAMBIT Overview](https://theiagen.notion.site/GAMBIT-7c1376b861d0486abfbc316480046bdc?pvs=4) |
        | Original Publication(s) | [GAMBIT (Genomic Approximation Method for Bacterial Identification and Tracking): A methodology to rapidly leverage whole genome sequencing of bacterial isolates for clinical identification](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0277575) <br> [TheiaEuk: a species-agnostic bioinformatics workflow for fungal genomic characterization](https://doi.org/10.3389/fpubh.2023.1198213) |