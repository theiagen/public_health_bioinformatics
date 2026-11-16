---
title: Task Fragment `medaka`
fragment: true
---
??? task "`Medaka`: Assembly Polishing (optional; default)"
    Medaka uses the raw reads to polish the assembly and generate a consensus sequence.

    !!! dna "Medaka supports legacy data models"
        If you have legacy data, the Medaka model will likely be `r941_min_hac_g507`.

    If the user wants to have Medaka auto-detect the model, set `medaka_auto_detect_model` to `true`. Alternatively, the user may specify a model with `medaka_model`. The default `medaka_model` is `r1041_e82_400bps_sup_v5.0.0`.

    The selected model is important for result accuracy. See also [Medaka's model documentation](https://github.com/nanoporetech/medaka?tab=readme-ov-file#models).

    When run in TheiaProk_ONT, the `--bacteria` flag is included; this behavior can be turned off with the `medaka_use_bacteria` optional input. In TheiaEuk_ONT, the `--bacteria` flag is not available.

    !!! techdetails "Medaka Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_medaka.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/polishing/task_medaka.wdl) |
        | Software Source Code | [Medaka on GitHub](https://github.com/nanoporetech/medaka) |
        | Software Documentation | [Medaka Documentation](https://github.com/nanoporetech/medaka#medaka) |
