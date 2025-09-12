<!-- if: assembly -->
??? task "`ShigEiFinder`: _Shigella_/EIEC Differentiation and Serotyping ==_using the assembly file as input_=="
<!-- endif -->
<!-- if: reads -->
??? task "`ShigEiFinder_reads`: _Shigella_/EIEC Differentiation and Serotyping ==_using Illumina read files as input_== (optional)"
    To activate the `shigeifinder_reads` task, set the `call_shigeifinder_reads_input` to be `true`. If set to `true`, `shigeifinder_reads` will run **in addition to** the assembly-based `shigeifinder` task.
<!-- endif -->

    ShigEiFinder differentiates _Shigella_ and enteroinvasive _E. coli_ (EIEC) using cluster-specific genes, identifies some serotypes based on the presence of O-antigen and H-antigen genes (_wzx_ and _wzy_), and predicts the number of virulence plasmids. It can serotype over 59 _Shigella_ and 22 EIEC serotypes using BLAST and BWA.

    !!! techdetails "ShigEiFinder Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_shigeifinder.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/escherichia_shigella/task_shigeifinder.wdl) |
        | Software Source Code | [ShigEiFinder on GitHub](https://github.com/LanLab/ShigEiFinder) |
        | Software Documentation | [ShigEiFinder on GitHub](https://github.com/LanLab/ShigEiFinder) |
        | Original Publication(s) | [Cluster-specific gene markers enhance Shigella and enteroinvasive Escherichia coli in silico serotyping](https://pubmed.ncbi.nlm.nih.gov/34889728/) |
