??? task "`StxTyper`: Identification and Typing of Shiga toxin (Stx) Genes ==_using the assembly file as input_=="
    StxTyper screens bacterial genome assemblies for shiga toxin genes and subtypes them into known subtypes and also looks for novel subtypes in cases where the detected sequences diverge from the reference sequences.
    
    Shiga toxin is the main virulence factor of Shiga-toxin-producing E. coli (STEC), though these genes are also found in Shigella species as well as some other genera more rarely, such as Klebsiella. [Please see this review paper that describes shiga toxins in great detail](https://doi.org/10.3390/microorganisms12040687).

    !!! tip "Running StxTyper via the TheiaProk workflows"
        The TheiaProk workflow will automatically run `stxtyper` on all E. coli and Shigella spp. samples, but _the user can opt-in to running the tool on any sample by setting the optional input variable `call_stxtyper` to `true` when configuring the workflow._
    
    Generally, `stxtyper` looks for _stxA_ and _stxB_ subunits that compose a complete operon. The A subunit is longer (in amino acid length) than the B subunit. Stxtyper attempts to detect these, compare them to a database of known sequences, and type them based on amino acid composition.  There typing algorithm and rules defining how to type these genes & operons will be described more completely in a publication that will be available in the future.
    
    The `stxtyper_report` output TSV is provided in [this output format.](https://github.com/ncbi/stxtyper/tree/v1.0.24?tab=readme-ov-file#output)

    This tool has been incorporated into v4.0.3 of AMRFinderPlus and runs behind-the-scenes when the user (or in this case, the TheiaProk workflow) provides the `amrfinder --organism Escherichia --plus` options.

    !!! techdetails "StxTyper Technical Details"

        |  | Links |
        | --- | --- |
        | Task | [task_stxtyper.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/escherichia_shigella/task_stxtyper.wdl) |
        | Software Source Code | [ncbi/stxtyper GitHub repository](https://github.com/ncbi/stxtyper) |
        | Software Documentation | [ncbi/stxtyper GitHub repository](https://github.com/ncbi/stxtyper) |
