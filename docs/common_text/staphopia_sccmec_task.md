??? task "`staphopia-sccmec`: Sequence Typing"
    SCC_mec_ (staphylococcal cassette chromosome _mec_) is a mobile genetic element of _Staphylococcus_ species. It includes a methicilin-resistant _mecA_ gene that is shared between _Staphylococcus_ strains via horizontal gene transfer, which leads to MRSA strains. SCC_mec_ has also been found to confer resistance to non-beta-lactam drugs as well, making it an important target for identifying antimicrobial resistance in _Staphylococcus aureus_.

    This tool assigns a SCC_mec_ type by using BLAST to compare the SCC_mec_ primers against the provided _S. aureus_ assembly. `staphopia-sccmec` reports `True` for exact primer matches and `False` all others. The [Hamming Distance](https://en.wikipedia.org/wiki/Hamming_distance) is also used in TheiaProk's implementation of the tool to include the number of mismatches so that the `False` results can be examined in more detail.

    !!! techdetails "staphopia-sccmec Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_staphopiasccmec.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/staphylococcus/task_staphopiasccmec.wdl) |
        | Software Source Code | [staphopia-sccmec on GitHub](https://github.com/staphopia/staphopia-sccmec) |
        | Software Documentation | [staphopia-sccmec on GitHub](https://github.com/staphopia/staphopia-sccmec) |
        | Original Publication(s) | [_Staphylococcus aureus_ viewed from the perspective of 40,000+ genomes](https://doi.org/10.7717/peerj.5261) |
