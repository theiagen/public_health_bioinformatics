??? task "`Flye`: _De novo_ Assembly"
    `flye_denovo` is a sub-workflow that performs _de novo_ assembly using Flye for ONT data and supports additional polishing and visualization steps.

    !!! tip "Ensure correct medaka model is selected if performing medaka polishing"
        In order to obtain the best results, the appropriate model must be set to match the sequencer's basecaller model; this string takes the format of {pore}\_{device}\_{caller variant}\_{caller_version}. See also <https://github.com/nanoporetech/medaka?tab=readme-ov-file#models>. If `flye` is being run on legacy data the medaka model will likely be `r941_min_hac_g507`. Recently generated data will likely be suited by the default model of `r1041_e82_400bps_sup_v5.0.0`.

    The detailed steps and tasks are as follows:

{{ include_md("common_text/porechop_task.md", indent=4, replacements={'??? task "`porechop`"' : '??? toggle "`Porechop`: Read Trimming (optional; off by default)"'}) }}

{{ include_md("common_text/flye_task.md", indent=4, replacements={'??? task "`flye`"' : '??? toggle "`Flye`: _De novo_ Assembly"'}) }}

{{ include_md("common_text/bandage_task.md", indent=4, replacements={"??? task": "??? toggle"}) }}

{{ include_md("common_text/polypolish_task.md", indent=4, replacements={"??? task": "??? toggle"}) }}

{{ include_md("common_text/medaka_task.md", indent=4, replacements={"??? task": "??? toggle"}) }}

{{ include_md("common_text/racon_task.md", indent=4, replacements={"??? task": "??? toggle"}) }}

{{ include_md("common_text/filter_contigs_task.md", indent=4, replacements={"??? task": "??? toggle"}, condition="flye") }}

{{ include_md("common_text/dnaapler_task.md", indent=4, replacements={"??? task": "??? toggle"}) }}

    !!! techdetails "Flye-Denovo Technical Details"
        |  | Links |
        | --- | --- |
        | Subworkflow | [wf_flye_denovo.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_flye_denovo.wdl) |
