??? task "`emmtyper`: Sequence Typing"
    The _Streptococcus pyogenes_ M protein (encoded by _emm_) is used for sequencing typing and disease surveillence as it is a major virulence factor. There are over 275 _emm_ types.

    `emmtyper` uses BLAST to compare the genome assembly to the CDC-curated trimmed _emm_ type database (by default). An _in silico_ PCR method is also available. The BLAST results are then processed to distinguish between _emm_ and _emm_-like alleles to derive the isolates' M-type. The predicted _emm_-type is reported in addition to any possible _emm_-like alleles and the functional _emm_ cluster.

    !!! techdetails "emm-typing-tool Technical Details"            
        |  | Links |
        | --- | --- |
        | Task | [task_emmtypingtool.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/streptococcus/task_emmtyper.wdl) |
        | Software Source Code | [emmtyper on GitHub](https://github.com/MDU-PHL/emmtyper/tree/master) |
        | Software Documentation | [emmtyper on GitHub](https://github.com/MDU-PHL/emmtyper/tree/master) |
        | Original Publication(s) | _emm typing scheme_: [A systematic and functional classification of Streptococcus pyogenes that serves as a new tool for molecular typing and vaccine development](https://doi.org/10.1093/infdis/jiu260) |
