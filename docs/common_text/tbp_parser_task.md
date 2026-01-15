??? task "`tbp-parser`: Interpretation and Parsing of TBProfiler JSON outputs (optional)"
    To activate this task, set `call_tbp_parser` to `true`

    [tbp-parser](https://github.com/theiagen/tbp-parser/) was developed by Theiagen in partnership with the California Department of Public Health. This tool adds useful drug resistance interpretation by applying expert rules and organizing the outputs from TBProfiler. To understand this module and its functions, [please examine the documentation for this tool](https://theiagen.github.io/tbp-parser/latest/). 

    This tool generates reports that can be automatically imported into your local LIMS system and is highly customizable with additional quality control metrics and fine-tuning. It is highly recommended to read the associated documentation to understand the full capabilities and configuration options available.

    _Please note that this tool has **not** been tested on ONT data and although it is available, result accuracy should be considered carefully._
    
    !!! techdetails "tbp-parser Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_tbp_parser.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/mycobacterium/task_tbp_parser.wdl) |
        | Software Source Code | [tbp-parser on GitHub](https://github.com/theiagen/tbp-parser/) |
        | Software Documentation | [tbp-parser Documentation](https://theiagen.github.io/tbp-parser) |
