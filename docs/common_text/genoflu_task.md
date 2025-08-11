??? task "`genoflu`"

    This task determines the whole-genome genotype of a H5N1 (currently only for the 2.3.4.5b clade of H5N1) flu sample by comparing each segment of the sample against a curated database of H5N1 references. Each segment is assigned a type, and the whole-genome genotype is assigned based on the combination of segment types, according to the [GenoFLU reference table](https://github.com/USDA-VS/GenoFLU/blob/main/docs/Genotyping_reference_for_US_H5_2.3.4.4b.pdf).
    
    !!! techdetails "GenoFLU Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_genoflu.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/orthomyxoviridae/task_genoflu.wdl) |
        | Software Source Code | [GenoFLU on GitHub](https://github.com/USDA-VS/GenoFLU) |
        | Software Documentation | [GenoFLU on GitHub](https://github.com/USDA-VS/GenoFLU) |
        | Original Publication(s) | [H5N1 highly pathogenic avian influenza clade 2.3.4.4b in wild and domestic birds: Introductions into the United States and reassortments, December 2021-April 2022](https://doi.org/10.1016/j.virol.2023.109860) |