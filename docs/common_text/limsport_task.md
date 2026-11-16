---
title: Task Fragment `limsport`
fragment: true
---
??? task "`limsport`: LIMS Report Generation and QC"

    LIMSport reads a downloaded Terra table to process the rows indicated in the `sample_names` variable using a configuration file. The configuration file is the heart of LIMSport.

    The [LIMSport tool documentation](https://theiagen.github.io/limsport/latest/) contains more details regarding the configuration file and other input parameters, and includes a config builder to help build your configuration file. [This builder is also included above for convenience](#config-builder).

    !!! techdetails "LIMSport Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_limsport.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/data_handling/task_limsport.wdl) |
        | Software Source Code | [LIMSport on GitHub](https://github.com/theiagen/limsport) |
        | Software Documentation | [LIMSport Documentation](https://theiagen.github.io/limsport/latest/) |
