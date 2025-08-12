??? task "`vadr`"

    VADR annotates and validates completed assembly files. For details on VADR default models/parameters, see the [*organism-specific parameters and logic* section](./theiacov.md#org-specific). VADR is also used to characterize individual Influenza segments and will extract each segment into its own fasta file.

    !!! techdetails "VADR Technical Details"

        |  | Links |
        | --- | --- |
        | Task | [task_vadr.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/advanced_metrics/task_vadr.wdl) |
        | Software Source Code | <https://github.com/ncbi/vadr> |
        | Software Documentation | <https://github.com/ncbi/vadr/wiki> |
        | Original Publication(s) | For SARS-CoV-2: *[Faster SARS-CoV-2 sequence validation and annotation for GenBank using VADR](https://doi.org/10.1093/nargab/lqad002)*<br> For non-SARS_CoV-2: [*VADR: validation and annotation of virus sequence submissions to GenBank*](https://doi.org/10.1186/s12859-020-3537-3) |