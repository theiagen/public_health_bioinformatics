??? task "`ShigaPass`: _Shigella_/EIEC differentiation and serotyping"
    ShigaPass predicts _Shigella_ spp. serotypes from assembled whole genomes, and can differentiate between _Shigella_, enteroinvasive _E. coli_ (EIEC), and non-_Shigella_/EIEC genomes.

    ShigaPass provides enhanced accuracy compared to other tools (e.g., ShigaTyper, and ShigeiFinder) and is recommended by the National Reference Center for _E. coli_, _Shigella_, and _Salmonella_ in France. 

    !!! warning "==ShigaPass works best with assemblies generated with SPAdes=="
        ShigaPass operates optimally using assemblies generated with SPAdes. Other assemblers, like SKESA (TheiaProk_Illumina's default), can lead to incorrect or incomplete results due to missing detetion of critical genes like _ipaH_. 

        Please be aware of this potential limitation in accuracy when interpreting your results.

    !!! tip "Assembler recommendations for ShigaPass"
        - **TheiaProk_Illumina_PE** & **TheiaProk_Illumina_SE**: We strongly recommend setting the `digger_denovo.assembler` input parameter to `"spades"` to ensure the best results from ShigaPass.
        - **TheiaProk_ONT**: By default, ONT assemblies are generated with `flye`. While these assemblies can be used successfully with ShigaPass, the results may be less accurate than if the data were assembled with SPAdes. Since SPAdes is not an option for ONT data, please be aware of this potential limitation in accuracy when interpreting your results.
        - **TheiaProk_FASTA**: Results may be less accurate if the data was not assembled with SPAdes. You may want to consider re-assembling your data with SPAdes for optimal performance of ShigaPass if possible. For publicly available data, you can check the original publication or SRA metadata to see if SPAdes was used for assembly.
    
    !!! techdetails "ShigaPass Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_shigapass.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/escherichia_shigella/task_shigapass.wdl) |
        | Software Source Code | [ShigaPass on GitHub](https://github.com/imanyass/ShigaPass) |
        | Software Documentation | [ShigaPass on GitHub](https://github.com/imanyass/ShigaPass) |
        | Original Publication(s) | [ShigaPass: an _in silico_ tool predicting _Shigella_ serotypes from whole-genome sequencing assemblies](https://doi.org/10.1099/mgen.0.000961) |
