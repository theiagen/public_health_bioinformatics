??? task "`Rasusa`: Read subsampling (optional, on by default)"

    The Rasusa task performs subsampling of the raw reads. By default, this task will subsample reads to a depth of 150X using the estimated genome length produced during the preceding raw read screen. The user can prevent the task from being launched by setting the `call_rasusa`variable to false. 

    The user can also provide an estimated genome length for the task to use for subsampling using the `genome_size` variable. In addition, the read depth can be modified using the `subsample_coverage` variable.
        
    !!! techdetails "Rasusa Technical Details"

        |  | Links |
        | --- | --- |
        | Task | [task_rasusa.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/task_rasusa.wdl) |
        | Software Source Code | [Rasusa on GitHub](https://github.com/mbhall88/rasusa) |
        | Software Documentation | [Rasusa on GitHub](https://github.com/mbhall88/rasusa) |
        | Original Publication(s) | [Rasusa: Randomly subsample sequencing reads to a specified coverage](https://doi.org/10.21105/joss.03941) |