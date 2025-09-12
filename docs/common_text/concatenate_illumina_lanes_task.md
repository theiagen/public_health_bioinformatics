??? task "`concatenate_illumina_lanes`: Concatenate Multi-Lane Illumina FASTQs"
    The `concatenate_illumina_lanes` task concatenates Illumina FASTQ files from multiple lanes into a single file. This task only runs if the `read1_lane2` input file has been provided. All read1 lanes are concatenated together and are used in subsequent tasks, as are the read2 lanes if applicable. These concatenated files are also provided as output.

    !!! techdetails "Concatenate Illumina Lanes Technical Details"
        The `concatenate_illumina_lanes` task is run before any downstream steps take place.
        
        |  | Links |
        | --- | --- |
        | Subworkflow | [wf_concatenate_illumina_lanes.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/file_handling/wf_concatenate_illumina_lanes.wdl)
