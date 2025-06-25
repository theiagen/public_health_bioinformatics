??? task "`nanoq`"

    Reads are filtered by length and quality using `nanoq`. By default, sequences with less than 500 basepairs and quality score lower than 10 are filtered out to improve assembly accuracy.

    !!! techdetails "Nanoq Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_nanoq.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_nanoq.wdl) |
        | Software Source Code | [Nanoq](https://github.com/esteinig/nanoq) |
        | Software Documentation | [Nanoq Documentation](https://github.com/esteinig/nanoq/blob/master/README.md) |
        | Original Publication(s) | [Nanoq: ultra-fast quality control for nanopore reads](https://doi.org/10.21105/joss.02991)