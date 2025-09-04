??? task "`SonneiTyper`: _Shigella sonnei_ identification, genotyping, and resistance mutation identification ==_for Illumina and ONT data only_=="
    SonneiTyper identifies _Shigella sonnei_, and uses **single-nucleotide variants** for genotyping and prediction of quinolone resistance in _gyrA_ (S83L, D87G, D87Y) and _parC_ (S80I). Outputs are provided in [a TSV format described here](https://github.com/katholt/sonneityping#example-output).

    SonneiTyper is a wrapper script around another tool, Mykrobe, that analyses the _S. sonnei_ genomes using the _S. sonnei_-specific genotyping scheme. SonneiTyper parses the Mykrobe predict results and tabulates the results.

    !!! techdetails "SonneiTyper Technical Details"

        |  | Links |
        | --- | --- |
        | Task | [task_sonneityping.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/escherichia_shigella/task_sonneityping.wdl) |
        | Software Source Code | [Mykrobe on GitHub](https://github.com/Mykrobe-tools/mykrobe)<br>[sonneityping on GitHub](https://github.com/katholt/sonneityping) |
        | Software Documentation | [Mykrobe Wiki](https://github.com/Mykrobe-tools/mykrobe/wiki)<br>[sonneityping on GitHub](https://github.com/katholt/sonneityping) |
        | Original Publication(s) | _Mykrobe tool_: [Antibiotic resistance prediction for Mycobacterium tuberculosis from genome sequence data with Mykrobe](https://doi.org/10.12688/wellcomeopenres.15603.1)<br>_S. sonnei genotyping scheme_: [Global population structure and genotyping framework for genomic surveillance of the major dysentery pathogen, _Shigella sonnei_](https://www.nature.com/articles/s41467-021-22700-4) |
