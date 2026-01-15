??? task "`vadr`"

    VADR (Viral Annotation DefineR) annotates and validates completed assembly files. For details on VADR default models/parameters, see the [*organism-specific parameters and logic* section](./theiacov.md#org-specific). It was primarily developed to test viral sequences to confirm they would be accepted to NCBI's GenBank data repository, but has found wide usage in general sequence validation and annotation.

    As part of the analysis of the assemblies, more than 70 types of unexpected characteristics, also known as _alerts_, can be reported. Any identified alerts can be found in the `vadr_alerts_list` output. Fatal alerts indicate that the sample is unlikely to be accepted to GenBank; non-fatal alerts are designated as passing sequences, but may still require further investigation. A full description of the potential alerts can be found on the [VADR README here](https://github.com/ncbi/vadr/blob/master/documentation/alerts.md), including details on how to allow sequencecs to pass despite having fatal alerts.

    !!! techdetails "VADR Technical Details"

        |  | Links |
        | --- | --- |
        | Task | [task_vadr.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/advanced_metrics/task_vadr.wdl) |
        | Software Source Code | <https://github.com/ncbi/vadr> |
        | Software Documentation | <https://github.com/ncbi/vadr/wiki> |
        | Original Publication(s) | For SARS-CoV-2: *[Faster SARS-CoV-2 sequence validation and annotation for GenBank using VADR](https://doi.org/10.1093/nargab/lqad002)*<br> For non-SARS_CoV-2: [*VADR: validation and annotation of virus sequence submissions to GenBank*](https://doi.org/10.1186/s12859-020-3537-3) |