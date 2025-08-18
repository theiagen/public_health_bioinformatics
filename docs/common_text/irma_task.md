<!-- if: assembly -->
??? toggle "`irma`: Assembly and Characterization ==_for flu in TheiaCoV_Illumina_PE & TheiaCoV_ONT_=="
<!-- endif -->
<!-- if: task -->
??? task "`irma`"
<!-- endif -->

    Cleaned reads are assembled using `irma` which stands for Iterative Refinement Meta-Assembler. IRMA first sorts reads to Flu genome segments using LABEL, then iteratively maps read to collection of reference sequences (in this case for Influenza virus) and iteratively edits the references to account for high population diversity and mutational rates that are characteristic of Influenza genomes. Assemblies produced by `irma` will be ordered from largest to smallest assembled flu segment. `irma` also performs typing and subtyping as part of the assembly process. Note: IRMA does not differentiate between Flu B Victoria and Yamagata lineages. For determining this information, please review the `abricate` task outputs which will provide this information.

<!-- if: assembly -->
    Due to the segmented nature of the Influenza genome and the various downstream bioinformatics tools that require the genome assembly, the IRMA task & TheiaCoV workflows output various genome assembly files. Briefly they are:

    - `assembly_fasta` - The full genome assembly in FASTA format, with 1 FASTA entry per genome segment. There should be 8 segments in total, but depending on the quality and depth of sequence data, some segments may not be assembled and nor present in this output file.
    - `irma_assembly_fasta_concatenated` - The full genome assembly in FASTA format, but with all segments concatenated into a single FASTA entry. This is not your typical FASTA file and is purposely created to be used with a custom Nextclade dataset for the H5N1 B3.13 genotype that is based on a concatenated reference genome.
    - `irma_<segment-abbreviation>_segment_fasta` - Individual FASTA files that only contain the sequence for 1 segment, for example the HA segment. There are 8 of these in total.

    General statistics about the assembly are generated with the `consensus_qc` task ([task_assembly_metrics.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_assembly_metrics.wdl)).
<!-- endif -->

    !!! techdetails "IRMA Technical Details" 
        |  | Links |
        | --- | --- |
        | Task | [task_irma.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_irma.wdl) |
        | Software Documentation | [IRMA website](https://wonder.cdc.gov/amd/flu/irma/) |
        | Original Publication(s) | [Viral deep sequencing needs an adaptive approach: IRMA, the iterative refinement meta-assembler](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-3030-6) |