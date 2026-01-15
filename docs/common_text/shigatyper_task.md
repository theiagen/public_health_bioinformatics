??? task "`ShigaTyper`: _Shigella_/EIEC differentiation and serotyping ==_for Illumina and ONT only_=="
    ShigaTyper predicts _Shigella_ spp. serotypes from Illumina or ONT read data. If the genome is not _Shigella_ or enteroinvasive _E. coli_ (EIEC), the results from this tool will state this. In the notes it provides, it also reports on the presence of _ipaB_  which is suggestive of the presence of the "virulent invasion plasmid".

    ShigaTyper works by mapping the sample sequence to a _Shigella_ reference sequence database using minimap2. Specifically, serotype prediction is made through the serotype-specific _wzx_ gene as O-antigen expression is dependent on this gene. Additional criteria are applied if the serotype could not be solely predicted by this gene; please see the publication for more details.
    
    !!! techdetails "ShigaTyper Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_shigatyper.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/escherichia_shigella/task_shigatyper.wdl) |
        | Software Source Code | [ShigaTyper on GitHub](https://github.com/CFSAN-Biostatistics/shigatyper) |
        | Software Documentation | [ShigaTyper on GitHub](https://github.com/CFSAN-Biostatistics/shigatyper) |
        | Original Publication(s) | [In Silico Serotyping Based on Whole-Genome Sequencing Improves the Accuracy of Shigella Identification](https://doi.org/10.1128/AEM.00165-19) |
