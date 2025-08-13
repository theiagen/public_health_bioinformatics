??? task "`flu_antiviral_substitutions`"

    This subworkflow determines if any antiviral mutations are present in the HA, NA, and MP segments of H1N1 or H3N2 flu sample, or any in non-subtype-specific PA, PB1, and PB2 segments.
            
    These mutations are identified by generating a multiple sequence alignment (MSA) between each individual flu segment and the respective reference genome using MAFFT. Amino acid mutations are then called from the MSA. The resulting mutations are compared against [a list of known amino-acid substitutions associated with antiviral resistance](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/drug_resistance/task_flu_antiviral_subs.wdl) and any matches are reported. 

    This list of amino-acid substitutions includes both substitutions reported in the scientific literature and those inferred to potentially cause antiviral resistance based on analogous antiviral mutations in other flu subtypes. [**A table with the explanation for each amino-acid substitution in the antiviral resistance task is available here**](../../assets/files/antiviral_resistance_flu_aa_substitutions_explanations.xlsx).

    The list of known amino-acid substitutions associated with resistance can be expanded via optional user input `antiviral_aa_subs` in the format "`NA:V95A,HA:I97V`", i.e. `Protein:AAPositionAA`. 

    ??? toggle "Currently, the default mutations considered confer resistance to the following antivirals"
        - A_315675
        - Amantadine
        - Compound_367
        - Favipiravir
        - Fludase
        - L_742_001
        - Laninamivir
        - Oseltamivir (tamiflu)
        - Peramivir
        - Pimodivir
        - Rimantadine
        - Xofluza
        - Zanamivir

    !!! techdetails "Antiviral Substitutions Technical Details"        
        |  | Links |
        | --- | --- |
        | Sub-workflow | [wf_influenza_antiviral_substitutions.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_influenza_antiviral_substitutions.wdl) |
        | Tasks | [task_mafft.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/alignment/task_mafft.wdl)<br>[task_flu_antiviral_subs.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/drug_resistance/task_flu_antiviral_subs.wdl) |
        | Original Publication(s) | [MAFFT Multiple Sequence Alignment Software Version 7: Improvements in Performance and Usability](https://doi.org/10.1093/molbev/mst010)<br>[Next-Generation Sequencing: An Eye-Opener for the Surveillance of Antiviral Resistance in Influenza](https://doi.org/10.1016/j.tibtech.2019.09.009)
