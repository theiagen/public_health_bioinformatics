
??? toggle "`Dnaapler`: Final Assembly Orientation"
    Dnaapler reorients contigs to start at specific reference points. Dnaapler supports the following modes, which can be indicated by filling the `dnaapler_mode` input variable with the desired mode. The default is `all`, which reorients contigs to start with `dnaA`, `terL`, `repA`, or `COG1474`.

    - **all**: Reorients contigs to start with `dnaA`, `terL`, `repA`, or `COG1474` (_Default_)
    - **chromosome**: Reorients to begin with the `dnaA` chromosomal replication initiator gene, commonly used for bacterial chromosome assemblies.
    - **plasmid**: Reorients to start with the `repA` plasmid replication initiation gene, ideal for plasmid assemblie
    - **phage**: Reorients to start with the `terL` large terminase subunit gene, used for bacteriophage assemblies
    - **archaea**: Reorients to start with the `COG1474` archaeal Orc1/cdc6 gene, relevant for archaeal assemblies
    - **custom**: Reorients based on a user-specified gene in amino acid FASTA format for experimental or unique workflows
    - **mystery**: Reorients to start with a random CDS for exploratory purposes
    - **largest**: Reorients to start with the largest CDS in the assembly, often useful for poorly annotated genomes
    - **nearest**: Reorients to start with the first CDS nearest to the sequence start, resolving CDS breakpoints
    - **bulk**: Processes multiple contigs to start with the desired start gene (`dnaA`, `terL`, `repA`, or custom)

    !!! techdetails "Dnaapler Technical Details"
        |  | Links |
        | --- | --- |
        | WDL Task | [task_dnaapler.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_dnaapler.wdl) |
        | Software Source Code | [Dnaapler on GitHub](https://github.com/gbouras13/dnaapler) |
        | Software Documentation | [Dnaapler Documentation](https://github.com/gbouras13/dnaapler?tab=readme-ov-file#dnaapler) |
        | Original Publication(s) | [Dnaapler: a tool to reorient circular microbial genomes](https://joss.theoj.org/papers/10.21105/joss.05968) |
