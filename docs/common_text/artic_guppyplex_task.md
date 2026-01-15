??? task "`artic_guppyplex`: Read Filtering"
    Reads are filtered by length with `artic_guppyplex`, which is a part of the [`ARTIC` protocol](https://artic.network/fieldbioinformatics/fieldbioinformatics-sop.html). Since TheiaCoV was developed primarily for amplicon-based viral sequencing, this task is included to remove chimeric reads that are either too short or too long.

    !!! techdetails "artic_guppyplex Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_artic_guppyplex.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_artic_guppyplex.wdl) |
        | Software Source Code | [ARTIC on GitHub](https://github.com/artic-network/fieldbioinformatics/) |
        | Software Documentation | [ARTIC Documentation](https://artic.readthedocs.io/en/latest/?badge=latest) |
