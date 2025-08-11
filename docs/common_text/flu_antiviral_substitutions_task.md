??? task "`flu_antiviral_substitutions`"

    This subworkflow and task determines which, if any, antiviral mutations are present in the sample. 
    
    The assembled HA, NA, PA, PB1 and PB2 segments are compared against [a list of known amino-acid substitutions associated with resistance](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/drug_resistance/task_flu_antiviral_subs.wdl) to the antivirals A_315675, Amantadine, compound_367, Favipiravir, Fludase, L_742_001, Laninamivir, Oseltamivir (tamiflu), Peramivir, Pimodivir, Rimantadine, Xofluza, and Zanamivir. The list of known amino-acid substitutions associated with resistance can be expanded via optional user input `antiviral_aa_subs` in the format "`NA:V95A,HA:I97V`", i.e. `Protein:AAPositionAA`. 

    The list of amino-acid substitutions associated with antiviral resistance includes both substitutions reported to confer antiviral resistance in the scientific literature and those inferred to potentially cause antiviral resistance based on an analogous mutation reported to confer antiviral resistance in another flu subtype. [A table with the explanation for each amino-acid substitution in the antiviral resistance task is available here](../../assets/files/antiviral_resistance_flu_aa_substitutions_explanations.xlsx).

    !!! techdetails "Antiviral Substitutions Technical Details"        
        |  | Links |
        | --- | --- |
        | Workflow | [wf_influenza_antiviral_substitutions.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_influenza_antiviral_substitutions.wdl) |
        | Task | [task_flu_antiviral_subs.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/drug_resistance/task_flu_antiviral_subs.wdl) |
        | Publication | [Next-Generation Sequencing: An Eye-Opener for the Surveillance of Antiviral Resistance in Influenza](https://doi.org/10.1016/j.tibtech.2019.09.009)
