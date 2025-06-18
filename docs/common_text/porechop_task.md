??? task "`porechop`"

<!-- if: theiaprok|theiaeuk  -->
    Read trimming is optional and can be enabled by setting the `run_porchop` input variable to true.
<!-- endif -->

    Porechop is a tool for finding and removing adapters from ONT data. Adapters on the ends of reads are trimmed, and when a read has an adapter in the middle, the read is split into two.

<!-- if: theiaviral -->
    The `porechop` task is optional and is turned off by default. It can be enabled by setting the `call_porechop` parameter to `true`.
<!-- endif -->

    !!! techdetails "Porechop Technical Details"
        |  | Links |
        | --- | --- |
        | WDL Task | [task_porechop.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_porechop.wdl) |
        | Software Source Code | [Porechop on GitHub](https://github.com/rrwick/Porechop) |
        | Software Documentation | [https://github.com/rrwick/Porechop#porechop](https://github.com/rrwick/Porechop#porechop) |