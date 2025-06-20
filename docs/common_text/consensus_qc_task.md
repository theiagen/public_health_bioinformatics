
??? task "`consensus_qc`"

    The consensus_qc task generates a summary of genomic statistics from a consensus genome. This includes the total number of bases, "N" bases, degenerate bases, and an estimate of the percent coverage to the reference genome.

    !!! techdetails "`consensus_qc` Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_consensus_qc.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_consensus_qc.wdl) |
        | Software Source Docker Image | [Theiagen Docker Builds: utility:1.1](https://github.com/theiagen/theiagen_docker_builds/blob/main/utility/1.1/Dockerfile) |