??? task "`checkv`"

    CheckV is a fully automated command-line pipeline for assessing the quality of viral genomes, including identification of host contamination for integrated proviruses, estimating completeness for genome fragments, and identification of closed genomes.

    By default, CheckV reports results on a contig-by-contig basis. The `checkv` task additionally reports both "weighted_contamination" and "weighted_completeness", which are average percents calculated across the total assembly that are weighted by contig length.

    !!! techdetails "CheckV Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_checkv.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/advanced_metrics/task_checkv.wdl) |
        | Software Source Code | [CheckV on Bitbucket](https://bitbucket.org/berkeleylab/checkv/src/master/) |
        | Software Documentation | [CheckV Documentation](https://bitbucket.org/berkeleylab/checkv/src/master/README.md) |
        | Original Publication(s) | [CheckV assesses the quality and completeness of metagenome-assembled viral genomes](https://doi.org/10.1038/s41587-020-00774-7) |
