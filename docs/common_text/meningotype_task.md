??? task "`meningotype`: _Neisseria meningitidis_ Serotyping"
    This tool performs _in silico_ typing of _N. meningitidis_. It performs the following functions: serogrouping, finetyping of _porA_ and _fetA_, _porB_ sequencing typing, and Bexsero Antigen Sequencing Typing (BAST) (_fHbp_, _NHBA_, _NadA_, and _PorA_). These results are all parsed and provided to the user as outputs.

    The allele database is extracted from PubMLST's _N. meningitidis_ database. This tool works by using BLAST to compare the sample sequence against the database entries.
    
    !!! techdetails "meningotype Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_meningotype.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/neisseria/task_meningotype.wdl) |
        | Software Source Code | [meningotype on GitHub](https://github.com/MDU-PHL/meningotype) |
        | Software Documentation | [meningotype on GitHub](https://github.com/MDU-PHL/meningotype) |
