??? task "`Racon`: Polishing of Flye assembly (alternative; optional)"
    Polishing is optional and can be skipped by setting the `skip_polishing` variable to true. If polishing is skipped, then neither Medaka or Racon will run.

    `Racon` is an alternative to using `medaka` for assembly polishing, and can be run by setting the `polisher` input to "racon".  Racon is a consensus algorithm designed for refining raw de novo DNA assemblies generated from long, uncorrected sequencing reads.

    !!! techdetails "Racon Technical Details"
        |  | Links |
        | --- | --- |
        | WDL Task | [task_racon.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/polishing/task_racon.wdl) |
        | Software Source Code | [Racon on GitHub](https://github.com/lbcb-sci/racon) |
        | Software Documentation | [Racon Documentation](https://github.com/lbcb-sci/racon#racon) |
        | Original Publication(s) | [Fast and accurate de novo genome assembly from long uncorrected reads](https://genome.cshlp.org/content/27/5/737) |
