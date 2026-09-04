---
title: Task Fragment `dorado`
fragment: true
---
??? task "`Dorado`: Assembly Polishing (alternative; optional)"
    Dorado uses raw sequencing reads to polish an assembly.

    !!! warning "Dorado DOES NOT support legacy data models"
        [Only the DNA R10.4.1 and RNA004 models are supported in Dorado](https://software-docs.nanoporetech.com/dorado/latest/models/list/). If you have legacy data generated with DNA R9.4.1 or RNA002, we recommend using Medaka instead.

    Dorado will by default auto-detect the basecalling model that was used to generate the data. If the user wants to specify their own model, they must set `dorado_auto_detect_model` to `false` and provide the model in the `dorado_model` input. Read groups will be ignored automatically; this behavior can be turned off with `dorado_ignore_read_groups`.

    The following models are available in the default Docker image in this task:

    ??? toggle "Available Dorado models"
          The following models are compatible with this version of Dorado:

          - dna_r10.4.1_e8.2_400bps_hac@v5.0.0_polish_rl
          - dna_r10.4.1_e8.2_400bps_hac@v5.0.0_polish_rl_mv
          - dna_r10.4.1_e8.2_400bps_sup@v5.0.0_polish_rl
          - dna_r10.4.1_e8.2_400bps_sup@v5.0.0_polish_rl_mv
          - dna_r10.4.1_e8.2_400bps_hac@v5.2.0_polish_rl
          - dna_r10.4.1_e8.2_400bps_hac@v5.2.0_polish_rl_mv
          - dna_r10.4.1_e8.2_400bps_sup@v5.2.0_polish_rl
          - dna_r10.4.1_e8.2_400bps_sup@v5.2.0_polish_rl_mv
          - dna_r10.4.1_e8.2_400bps_hac@v6.0.0_polish_rl
          - dna_r10.4.1_e8.2_400bps_hac@v6.0.0_polish_rl_mv
          - dna_r10.4.1_e8.2_400bps_polish_bacterial_methylation_v5.0.0
          - dna_r10.4.1_e8.2_400bps_hac@v4.2.0_polish
          - dna_r10.4.1_e8.2_400bps_sup@v4.2.0_polish
          - dna_r10.4.1_e8.2_400bps_hac@v4.3.0_polish
          - dna_r10.4.1_e8.2_400bps_sup@v4.3.0_polish

          Specify a model at your own risk.

    The selected model is important for result accuracy. See also [Dorado's model documentation](https://software-docs.nanoporetech.com/dorado/latest/models/models/).

    When run in TheiaProk_ONT, the `--bacteria` flag is included; this behavior can be turned off with the `dorado_use_bacteria` optional input. In TheiaEuk_ONT, the `--bacteria` flag is not available.

    !!! techdetails "Dorado Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_dorado.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/polishing/task_dorado.wdl) |
        | Software Source Code | [Dorado on GitHub](https://github.com/nanoporetech/dorado) |
        | Software Documentation | [Dorado Documentation](https://software-docs.nanoporetech.com/dorado/latest/) |
