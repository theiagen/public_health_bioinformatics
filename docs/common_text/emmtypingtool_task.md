??? task "`emm-typing-tool`: Sequence Typing ==_for Illumina_PE only_=="
    The _Streptococcus pyogenes_ M protein (encoded by _emm_) is used for sequencing typing and disease surveillence as it is a major virulence factor. There are over 275 _emm_ types.

    `emm-typing-tool` maps the reads to a CDC-curated _emm_ type database using bowtie2 to identify any _emm_ genes. Alleles with 100% coverage and over 90% identity are selected, and the allele with the highest percent identity is generally reported (see the [decision tree](https://github.com/ukhsa-collaboration/emm-typing-tool/blob/master/decision_algorithm.png) for the nuances).

    !!! techdetails "emm-typing-tool Technical Details"            
        |  | Links |
        | --- | --- |
        | Task | [task_emmtypingtool.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/streptococcus/task_emmtypingtool.wdl) |
        | Software Source Code | [emm-typing-tool on GitHub](https://github.com/ukhsa-collaboration/emm-typing-tool) |
        | Software Documentation | [emm-typing-tool on GitHub](https://github.com/ukhsa-collaboration/emm-typing-tool) |
        | Original Publication(s) | _emm typing scheme_: [A systematic and functional classification of Streptococcus pyogenes that serves as a new tool for molecular typing and vaccine development](https://doi.org/10.1093/infdis/jiu260) |
