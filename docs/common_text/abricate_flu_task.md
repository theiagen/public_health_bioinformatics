??? task "`abricate`"

    ABRicate assigns types and subtype/lineages for flu samples using a version of the [INSaFLU](https://github.com/INSaFLU/INSaFLU) ("INSide the FLU") database [described here](https://github.com/epi2me-labs/wf-flu/tree/master/data/primer_schemes/V1/blastdb/insaflu). 

    ABRicate typically works by screening contigs for the presence of acquired resistance genes, but when using the INSaFLU database, the algorithm works by assigning contigs to the most closely corresponding viral segment in the INSaFLU database, which is used to call the flu type and subtype.
    
    !!! techdetails "ABRicate Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_abricate.wdl (abricate_flu subtask)](https://github.com/theiagen/public_health_bioinformatics/blob/2dff853defc6ea540a058873f6fe6a78cc2350c7/tasks/gene_typing/drug_resistance/task_abricate.wdl#L59) |
        | Software Source Code | [ABRicate on GitHub](https://github.com/tseemann/abricate) |
        | Software Documentation | [ABRicate on GitHub](https://github.com/tseemann/abricate) |
        | Original Publication(s) | _INSaFLU database_: [INSaFLU: an automated open web-based bioinformatics suite "from-reads" for influenza whole-genome-sequencing-based surveillance](https://doi.org/10.1186/s13073-018-0555-0) |
