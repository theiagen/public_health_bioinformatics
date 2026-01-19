??? task "`bbmap_reformat`: FASTQ Reformatting"
    The `bbmap_reformat` task takes in interleaved FASTQ files, repairs them if necessary, and subsequently deinterleaves the repaired or input FASTQ files.

<!-- if: theiacov -->
    This task is run when `irma_return_aligned_reads = true` and acts on the concatenated segment FASTQ files returned from IRMA. 
<!-- endif -->

    !!! techdetails "`bbmap_reformat` Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_bbmap_reformat.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/data_handling/task_bbmap_reformat.wdl) |
        | Software Source Code | [BBMap on SourceForge](https://sourceforge.net/projects/bbmap/) |
        | Software Documentation | [BBDuk Guide (archived)](https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) |