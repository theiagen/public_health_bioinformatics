---
title: Task Fragment `dorado_trim`
fragment: true
---
??? task "`dorado_trim`: Custom Primer Trimming (optional)"
    If the optional input `custom_primers` is provided, this task is activated to trim the primer sequences from the beginning and end of the demultiplexed reads.

    To determine how to format the FASTA file that is expected in `custom_primers` please see the [Dorado documentation](https://dorado-docs.readthedocs.io/en/latest/barcoding/custom_primers/), specifically the section on "Custom adapter/primer file format".

    !!! tip "Older Dorado Version Used"
        The Dorado version used in this task is not the most up-to-date version (set to v0.8.3) due to a bug in the Dorado subcommand in the latest version (v0.9.0). This will be updated in the future when the bug has been resolved by the Dorado developers.

    !!! techdetails "Dorado Trimming Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_dorado_trim.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/basecalling/task_dorado_trim.wdl) |
        | Software Source Code | [Dorado on GitHub](https://github.com/nanoporetech/dorado/) |
        | Software Documentation | [Dorado ReadTheDocs](https://dorado-docs.readthedocs.io/en/latest/) |
