# Glossary

This glossary defines the acronyms, file formats, and technical terms used throughout the PHB documentation. It is written for a public health audience and favors plain-language explanations over strict formal definitions.

!!! dna "Hover for definitions in the docs"
    Organization abbreviations (e.g. `ENA`, `GISAID`, `INSDC`) appear as **hover tooltips** throughout the documentation — hover over one (or focus it with the keyboard) to see a short definition without leaving the page. This page is the complete reference for every term below.

    Contributors: to add or edit an entry, see the [Documentation Contribution Guide](../contributing/doc_contribution.md#glossary-and-abbreviations).

## Sequencing and analysis terms

adapter
:   A short synthetic DNA sequence added to fragments during library preparation; adapters are removed from reads before analysis (see read trimming).

allele
:   One of the alternative versions of a gene or genetic marker.

amplicon
:   A piece of DNA produced by amplifying a specific target region, for example by PCR.

ANI
:   Average nucleotide identity — a percentage measure of how similar two genomes are, often used to confirm species identity.

assembly
:   The reconstructed genome sequence produced by piecing together sequencing reads.

barcode
:   A short known sequence added to each sample so that pooled samples can be told apart (demultiplexed) after sequencing.

base pair (bp)
:   A single rung of the DNA double helix; the standard unit for measuring sequence length (e.g., a 5,000 bp genome).

basecalling
:   Converting the raw signal from a sequencer (especially Oxford Nanopore) into the actual sequence of bases (A, C, G, T).

CDS
:   Coding sequence — the part of a gene that is translated into protein.

consensus sequence
:   A single representative sequence derived from many overlapping reads, taking the most common base at each position.

contig
:   A contiguous stretch of DNA sequence assembled from overlapping reads.

coverage (depth)
:   How many times, on average, each base in a genome is represented by sequencing reads. Higher coverage generally means more reliable results.

de novo assembly
:   Assembling a genome from reads alone, without aligning to a reference genome. (_De novo_ is Latin for "from the beginning".)

genome annotation
:   Identifying and labeling the features (such as genes and coding sequences) within an assembled genome.

in silico
:   Performed by computer simulation rather than in the physical laboratory (for example, "in silico PCR").

indel
:   An insertion or deletion of one or more bases in a DNA sequence.

k-mer
:   A short subsequence of DNA of a fixed length _k_, used for fast sequence comparison and classification.

locus (plural: loci)
:   A fixed position in a genome, such as the location of a particular gene.

MAGs
:   Metagenome-assembled genomes — genomes reconstructed from mixed (metagenomic) sequencing data.

metagenomics
:   Sequencing all of the DNA in a mixed sample (containing many organisms) rather than a single isolate.

MNPs
:   Multi-nucleotide polymorphisms — variants affecting several adjacent bases at once.

MSA
:   Multiple sequence alignment — aligning three or more sequences to compare shared positions.

N50
:   An assembly quality metric: the contig length such that half of the assembly sits in contigs of that length or longer. Higher is generally better.

pangenome
:   The full set of genes found across a group of related genomes — the core genome (genes shared by all) plus the accessory genome (genes present in only some).

PCR
:   Polymerase chain reaction — a laboratory method for amplifying (copying) DNA.

PE / SE
:   Paired-end (PE) sequencing reads are generated from both ends of a DNA fragment; single-end (SE) reads are generated from one end only.

phylogenetic tree
:   A branching diagram that shows the inferred evolutionary relationships among a set of sequences.

quality score (Phred score)
:   A measure of how confident the sequencer is in each called base; higher scores mean fewer expected errors.

read
:   A single sequence of bases produced by a sequencing instrument.

read mapping
:   Aligning sequencing reads to a reference genome to determine where each read originated.

read trimming
:   Removing low-quality bases and adapter sequences from reads before assembly or analysis.

recombination
:   The exchange of genetic material between organisms; recombinant regions are often removed before building bacterial phylogenies (e.g., by Gubbins).

reference genome
:   A previously assembled, high-quality genome used as a standard for comparison when aligning reads or calling variants.

SNP / SNV
:   A single nucleotide polymorphism (SNP) or single nucleotide variant (SNV) is a single-letter difference in DNA between samples.

tNGS
:   Targeted next-generation sequencing — sequencing focused on specific genomic regions of interest.

variant calling
:   Identifying the positions where a sample's sequence differs from a reference genome (such as SNPs and indels).

WGS
:   Whole-genome sequencing — sequencing the entire genome of an organism.

## File formats

BAM / SAM
:   Sequence Alignment Map (SAM) is a text file of sequencing reads aligned to a reference genome; BAM is its compressed binary equivalent.

BED
:   Browser Extensible Data — a simple file format that lists genomic regions by their start and end coordinates.

CSV / TSV
:   Comma-separated values (CSV) and tab-separated values (TSV) are plain-text table formats where each line is a row and columns are separated by commas or tabs.

FASTA
:   A text format for storing nucleotide or protein sequences.

FASTQ
:   A text format for storing sequencing reads together with a quality score for each base.

GBFF
:   GenBank Flat File — an annotated sequence file format from NCBI, combining a genome sequence with its feature annotations.

GenBank
:   An NCBI sequence format that pairs a nucleotide sequence with its annotations; the flat-file version is the GBFF.

GFF / GFF3
:   General Feature Format — a file describing features (genes, coding regions, etc.) within a genome. GFF3 is version 3 of the format.

JSON
:   JavaScript Object Notation — a structured text format for storing data; used, for example, by Auspice to display phylogenetic trees.

Newick
:   A compact text format for representing phylogenetic trees (files often end in `.nwk`).

POD5
:   The raw signal file format produced by Oxford Nanopore sequencers.

VCF
:   Variant Call Format — a file that lists the genetic variants (such as SNPs and indels) found in a sample relative to a reference.

## Genomic characterization terms

AMR
:   Antimicrobial resistance — the ability of a microbe to survive drugs that would normally kill it or stop its growth.

antigen
:   A molecule (often on a microbe's surface) that the immune system recognizes; surface antigens are the basis of serotyping (e.g., the O and H antigens of _E. coli_).

ARG
:   Antimicrobial resistance gene — a gene that confers resistance to one or more antimicrobial drugs.

biotype
:   A subdivision of a species distinguished by physiological or biochemical traits (for example, the "Classical" and "El Tor" biotypes of _Vibrio cholerae_).

clade
:   A group of organisms that all descend from a single common ancestor on a phylogenetic tree.

genotype
:   The genetic makeup of an organism, or a group defined by shared genetic features.

GPSC
:   Global Pneumococcal Sequence Cluster — a standardized scheme for defining and naming pneumococcal strains.

lineage
:   A group of organisms connected by unbroken descent from a common ancestor; for viruses like SARS-CoV-2, lineages are named (for example, by Pangolin).

MLST / cgMLST / wgMLST
:   Multi-locus sequence typing classifies strains by the alleles present at several genes. Core genome (cgMLST) and whole genome (wgMLST) MLST extend this to hundreds or thousands of genes.

PBP
:   Penicillin-binding protein — bacterial proteins targeted by beta-lactam antibiotics; used for pneumococcal typing.

plasmid
:   A small, circular piece of DNA separate from the chromosome that can move between bacteria and often carries resistance or virulence genes.

point mutation
:   A change at a single position in the DNA; some antimicrobial resistance is caused by point mutations rather than acquired genes.

serogroup
:   A group of related serotypes (for example, the O1 and O139 serogroups of _Vibrio cholerae_).

serotype
:   A subdivision of a species distinguished by molecules on its surface (antigens).

serovar
:   A serologically distinct variant within a species; used interchangeably with serotype, especially for _Salmonella_.

ST
:   Sequence type — a strain designation assigned by MLST.

strain
:   A specific isolate or genetically distinct variant within a species.

subtype
:   A subdivision within a species or type; for influenza, subtypes are defined by the HA and NA surface proteins (e.g., H1N1, H3N2).

taxon (plural: taxa)
:   A named group of organisms, such as a species or genus.

virulence
:   The capacity of a microbe to cause disease; virulence genes contribute to this capacity.

XDR
:   Extensively drug-resistant — resistant to nearly all available antimicrobials.

## Platforms and tools

Auspice
:   The web-based visualization tool of Nextstrain, used to interactively explore phylogenetic trees — for example, the JSON output of the Augur workflow — at [auspice.us](https://auspice.us).

Nextstrain
:   An open-source project and toolkit for real-time tracking of pathogen evolution; its Augur (analysis) and Auspice (visualization) components underlie many phylogenetic builds.

ONT
:   Oxford Nanopore Technologies — a long-read sequencing platform.

## Databases, repositories, and organizations

ARLN
:   Antimicrobial Resistance Laboratory Network — a CDC network of public health laboratories for AMR surveillance.

BioProject
:   An NCBI record that groups the data from a single study or research effort under one accession.

BioSample
:   An NCBI record describing the biological source material for a sequenced sample; each is assigned a BioSample accession.

CDC
:   Centers for Disease Control and Prevention — the U.S. national public health agency.

DDBJ
:   DNA Data Bank of Japan — the Japanese member of the INSDC and source of the DRA.

GISAID
:   Global Initiative on Sharing All Influenza Data — a public repository for viral genomic data.

GTDB
:   Genome Taxonomy Database — a standardized microbial taxonomy based on genome data.

INSDC
:   International Nucleotide Sequence Database Collaboration — the partnership between NCBI (SRA), ENA, and DDBJ (DRA) that shares sequence data globally.

NCBI
:   National Center for Biotechnology Information — a major U.S. host of public biological databases.

PHA4GE
:   Public Health Alliance for Genomic Epidemiology — an international group that develops public health genomics standards and best practices.

RefSeq
:   NCBI's curated Reference Sequence database of well-annotated genomes.

SRA / ENA / DRA
:   The Sequence Read Archive (SRA, at NCBI), European Nucleotide Archive (ENA), and DDBJ Sequence Read Archive (DRA) are the three public repositories of raw sequencing data that make up the INSDC.

## Pathogens and organisms

EIEC
:   Enteroinvasive _Escherichia coli_ — a pathogenic _E. coli_ group closely related to, and differentiated from, _Shigella_.

HIV
:   Human immunodeficiency virus.

MPXV
:   Monkeypox (mpox) virus.

MRSA
:   Methicillin-resistant _Staphylococcus aureus_.

MTBC
:   _Mycobacterium tuberculosis_ complex — the group of closely related species that cause tuberculosis.

RSV
:   Respiratory syncytial virus.

SARS-CoV-2
:   Severe acute respiratory syndrome coronavirus 2 — the virus that causes COVID-19.

WNV
:   West Nile virus.

## Acronyms and abbreviations

CLI
:   Command-line interface — running software by typing commands rather than using a graphical interface.

HRRT
:   Human read removal tool — software that removes human sequences from a dataset before sharing.

PHB
:   Public Health Bioinformatics — this collection of bioinformatics workflows.

QC
:   Quality control — checking that data meet defined quality thresholds before use.

SOP
:   Standard operating procedure — a step-by-step protocol for performing a task consistently.

WDL
:   Workflow Description Language — the language used to write the PHB workflows.
