??? task "`hicap`: Serotyping"
    `hicap` identifies the _cap_ locus serotype in _Haemophilus influenzae_ assemblies. As described in the `hicap` documentation:

    > The _cap_ locus of _H. influenzae_ is categorised into 6 different groups based on serology (a-f). There are three functionally distinct regions of the _cap_ locus, designated `region I`, `region II`, and `region III`. Genes within `region I` (`bexABCD`) and `region III` (`hcsAB`) are associated with transport and post-translation modification. The `region II` genes encode serotype-specific proteins, with each serotype (a-f) having a distinct set of genes. _cap_ loci are often subject to structural changes (e.g. duplication, deletion) making the process of _in silico_ typing and characterisation of loci difficult.

    !!! techdetails "hicap Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_hicap.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/haemophilus/task_hicap.wdl) |
        | Software Source Code | [hicap on GitHub](https://github.com/scwatts/hicap) |
        | Software Documentation | [hicap on GitHub](https://github.com/scwatts/hicap) |
        | Original Publication(s) | [hicap: In Silico Serotyping of the Haemophilus influenzae Capsule Locus](https://doi.org/10.1128/jcm.00190-19) |
