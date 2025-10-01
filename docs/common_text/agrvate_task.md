??? task "`AgrVATE`: Sequence Typing"
    This tool identifies the accessory gene regulator (_agr_) locus type and reports possible variants in the _agr_ operon. The _agr_ system is a major regulator of virulence phenotypes in _Staphylococcus aureus_. Many _S. aureus_ strains often have nonfunctional _agr_ activity due to various loss-of-function mutations in the _agr_ operon. This may be associated with increased disease severity, making its detection clinically significant.

    AgrVATE workings by using BLAST to compare a database of unique _agr_-group-specific k-mers against the provided assembly. The _agr_ operon is then extracted using _in-silico_ PCR computational methods. If an _agr_-group is assigned, variants are called using an Agr-group specific reference operon in order to detect possible non-functional _agr_. If the _agr_-group was untypeable, the sequence is searched for potential non-canonical _agrD_.

    !!! techdetails "AgrVATE Technical Details"

        |  | Links |
        | --- | --- |
        | Task | [task_agrvate.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/staphylococcus/task_agrvate.wdl) |
        | Software Source Code | [AgrVATE on GitHub](https://github.com/VishnuRaghuram94/AgrVATE) |
        | Software Documentation | [AgrVATE on GitHub](https://github.com/VishnuRaghuram94/AgrVATE) |
        | Original Publication(s) | [Species-Wide Phylogenomics of the _Staphylococcus aureus_ Agr Operon Revealed Convergent Evolution of Frameshift Mutations](https://doi.org/10.1128/spectrum.01334-21) |
