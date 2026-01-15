??? task "`pasty`: Serotyping"
    `pasty` is a tool for _in silico_ serogrouping of _Pseudomonas aeruginosa_ isolates. `pasty` was developed by Robert Petit, based on the [PAst](https://github.com/Sandramses/PAst) tool from the Centre for Genomic Epidemiology.

    `pasty` uses the [`camlhmp`](https://github.com/rpetit3/camlhmp) tool to identify the serogroup. It uses BLAST to compare an input asssembly against a set of O-antigens. The serogroup can be predicted based off of those results.

    !!! techdetails "pasty Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_pasty.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/pseudomonas/task_pasty.wdl) |
        | Software Source Code | [pasty on GitHub](https://github.com/rpetit3/pasty) |
        | Software Documentation | [pasty on GitHub](https://github.com/rpetit3/pasty) |
        | Original Publication(s) | [Application of Whole-Genome Sequencing Data for O-Specific Antigen Analysis and In Silico Serotyping of Pseudomonas aeruginosa Isolates.](https://journals.asm.org/doi/10.1128/JCM.00349-16) |
