??? task "`abricate_vibrio`: _Vibrio_ Characterization"
    The `abricate_vibrio` task is used to perform general characterization of _Vibrio_ genomes using a database of target sequences that are traditionally used in PCR methods. The sequences included in the database are as follows:

    ??? toggle "Resistence Gene Database"
        | Sequence Name | Sequence Role | Purpose in Database |
        | --- | --- | --- |
        | _toxR_ | Transcriptional activator | Species marker where presence identifies _V. cholerae_  |
        | _ompW_ | Outer Membrane Protein | Species marker where presence identifies _V. cholerae_  |
        | _ctxA_ | Cholera toxin | Indicates cholera toxin production |
        | _tcpA_classical_ | Toxin co-pilus A allele, associated with the Classical biotype | Used to infer identity as Classical biotype |
        | _tcpA_ElTor_ | Toxin co-pilus A allele, associated with the El Tor biotype | Used to infer identity as El Tor biotype |
        | _wbeN_ | O antigen encoding region | Used to infer identity as O1 serogroup |
        | _wbfR_ | O antigen encoding region | Used to infer identity as O139 serogroup |

        This database was developed via communication with Dr. Christine Lee, of the National Listeria, Yersinia, Vibrio and Enterobacterales Reference Laboratory within the Enteric Diseases Laboratory Branch at CDC. It is identical to the database used in the `srst2` task except it is formatted for ABRicate.

    This task works by using the ABRicate tool, which uses BLAST to compare the assembled genome to the target sequences in the database and then reporting the details of the genes that pass quality thresholds. The presence or absence of specific genes are used to verify the species, identify cholera toxin production, and designation both the biotype and serogroup of the sample. See the table above for the genes used for each of these purposes.

    !!! techdetails "abricate_vibrio Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_abricate.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/vibrio/task_abricate_vibrio.wdl) |
        | Software Source Code | [ABRicate on GitHub](https://github.com/tseemann/abricate) |
        | Software Documentation | [ABRicate on GitHub](https://github.com/tseemann/abricate) |
