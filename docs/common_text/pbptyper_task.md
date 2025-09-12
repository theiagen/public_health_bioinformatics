??? task "`PBPtyper`: Penicillin-Binding Protein Genotyping"
    The Penicillin-binding proteins (PBP) are responsible for the minimum inhibitory concentration (MIC) phenotypes for beta-lactam antibiotics, which are the drugs of choice for treating pneumococcal infections. In _Streptococcus pneumoniae_, these PBP genes (PBP1a, PBP2b, and PBP2x) can be identified and typed with PBPTyper to help predict beta-lactam resistance levels in _S. pneumoniae_, which can be invaluable in disease surveillance.

    pbptyper uses BLAST and average nucleotide identity (ANI) to identify the closest matching PBP types from a database of known PBP types and uses the PBP typing scheme described by the linked publication below.
    
    !!! techdetails "pbptyper Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_pbptyper.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/streptococcus/task_pbptyper.wdl) |
        | Software Source Code | [pbptyper on GitHub](https://github.com/rpetit3/pbptyper) |
        | Software Documentation | [pbptyper on GitHub](https://github.com/rpetit3/pbptyper) |
        | Original Publication(s) | _PBP typing method_: [Penicillin-Binding Protein Transpeptidase Signatures for Tracking and Predicting β-Lactam Resistance Levels in _Streptococcus pneumoniae_](https://journals.asm.org/doi/full/10.1128/mBio.00756-16) |
