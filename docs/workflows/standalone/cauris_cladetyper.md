# Cauris_CladeTyper

!!! warning "NEEDS WORK!!!!"
    This page is under construction and will be updated soon.

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Standalone](../../workflows_overview/workflows_type.md/#standalone) | [Mycotics](../../workflows_overview/workflows_kingdom.md#mycotics) | PHB v1.0.0 | Yes | Sample-level |

## Cauris_CladeTyper_PHB

The Cauris_CladeTyper_PHB Workflow is designed to assign clade to _Candida auris_ Whole Genome Sequencing assemblies based on their genomic sequence similarity to the five clade-specific reference files. Clade typing is essential for understanding the epidemiology and evolutionary dynamics of this emerging multidrug-resistant fungal pathogen.

### Inputs

### Workflow Tasks

The Cauris_Cladetyper Workflow for _Candida auris_ employs GAMBIT for taxonomic identification, comparing whole genome sequencing data against reference databases to accurately classify _Candida auris_ isolates. A custom database featuring five clade-specific _Candida auris_ reference genomes facilitates clade typing. Sequences undergo genomic signature comparison against the custom database, enabling assignment to one of the five _Candida auris_ clades (Clade I to Clade V) based on sequence similarity and phylogenetic relationships. This integrated approach ensures precise clade assignments, crucial for understanding the genetic diversity and epidemiology of _Candida auris_.

### Outputs
