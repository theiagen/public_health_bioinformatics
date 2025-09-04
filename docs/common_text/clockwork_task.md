??? task "`Clockwork`: Read Decontamination ==_for Illumina PE only_=="
    [Clockwork](https://github.com/iqbal-lab-org/clockwork/wiki) decontaminates paired-end _Mycobacterium tuberculosis_ data by removing all non-TB reads. At a high level, the sample is processed by aligning the reads with BWA to the H37Rv reference genome, and retaining only reads that have been mapped. This greatly improves the quality of any called variants and ensures that any variants called by TBProfiler are of suitable reliability.

    !!! techdetails "Clockwork Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_clockwork.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/mycobacterium/task_clockwork.wdl) |
        | Software Source Code | [Clockwork on GitHub](https://github.com/iqbal-lab-org/clockwork) |
        | Software Documentation | [Clockwork Wiki](https://github.com/iqbal-lab-org/clockwork/wiki) |
        | Original Publication(s) | _Clockwork tool_: [Minos: variant adjudication and joint genotyping of cohorts of bacterial genome](https://doi.org/10.1186/s13059-022-02714-x) |
