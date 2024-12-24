# Cauris_CladeTyper

!!! warning "NEEDS WORK!!!!"
    This page is under construction and will be updated soon.

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Standalone](../../workflows_overview/workflows_type.md/#standalone) | [Mycotics](../../workflows_overview/workflows_kingdom.md#mycotics) | PHB vX.X.X | Yes | Sample-level |

## Cauris_CladeTyper_PHB

The Cauris_CladeTyper_PHB Workflow is designed to assign the clade to _Candida auris_ (also known as _Candidozyma auris_) WGS assemblies based on their genomic sequence similarity to the five clade-specific reference files. Clade typing is essential for understanding the epidemiology and evolutionary dynamics of this emerging multidrug-resistant fungal pathogen.

### Inputs

### Workflow Tasks

!!! task "Cauris_Cladetyper"
    The Cauris_Cladetyper Workflow for _Candida auris_ employs GAMBIT for taxonomic identification, comparing whole genome sequencing data against reference databases to accurately classify _Candida auris_ isolates.

    A custom GAMBIT database is created using five clade-specific _Candida auris_ reference genomes. Sequences undergo genomic signature comparison against this database, which then enables assignment to one of the five _Candida auris_ clades (Clade I to Clade V) based on sequence similarity and phylogenetic relationships. This integrated approach ensures precise clade assignments, crucial for understanding the genetic diversity and epidemiology of _Candida auris_.

    See more information on the reference information for the five clades below:

    | Clade | Assembly/RefSeq Accession | Assembly Name | Strain | BioSample Accession
    |---|---|---|---|
    | Clade I | GCA_002759435.2 | Cand_auris_B8441_V2 | B8441 | SAMN05379624 |
    | Clade II | GCA_003013715.2 | ASM301371v2 | B11220 | SAMN05379608 |
    | Clade III | GCA_002775015.1 | Cand_auris_B11221_V1 | B11221 | SAMN05379609 |
    | Clade IV | GCA_003014415.1 | Cand_auris_B11243 | B11243 | SAMN05379619 |
    | Clade V | GCA_016809505.1 | ASM1680950v1 | | SAMN11570381 |

### Outputs


