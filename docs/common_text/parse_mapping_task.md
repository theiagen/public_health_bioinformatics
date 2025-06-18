??? task "`parse_mapping`"

<!-- if: sam_to_sorted_bam-->
    The `sam_to_sorted_bam` sub-task converts the output SAM file from the `minimap2` task and converts it to a BAM file. It then sorts the BAM file by coordinate, and creates a BAM index file. This processed BAM is required for the `clair3` variant calling task.

    ??? dna "`min_map_quality`"
        This parameter accepts an integer value to set the minimum mapping quality for variant calling and subsequent consensus sequence generation. The default value is `20`.
<!-- endif -->

<!-- if: bam_to_unaligned_fastq-->
    The `bam_to_unaligned_fastq` task will extract a FASTQ file of reads that failed to align, while removing unpaired reads. 
<!-- endif -->

<!-- if: theiaviral_mask_low_coverage-->
    The `mask_low_coverage` sub-task is used to mask low coverage regions in the `reference_fasta` file to improve the accuracy of the final consensus genome. Coverage thresholds are defined by the `min_depth` parameter, which specifies the minimum read depth required for a base to be retained. Bases falling below this threshold are replaced with "N"s to clearly mark low confidence regions. The masked reference is then combined with variants from the `clair3` task to produce the final consensus genome.

    ??? dna "`min_depth`"
        This parameter accepts an integer value to set the minimum read depth for variant calling and subsequent consensus sequence generation. The default value is `10`.
<!-- endif -->

    !!! techdetails "`parse_mapping` Technical Details"
        | | Links |
        |---|---|
        | Task | [task_parse_mapping.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/data_handling/task_parse_mapping.wdl) |
        | Software Source Code | [samtools on GitHub](https://github.com/samtools/samtools) |
        | Software Documentation | [samtools](https://www.htslib.org/doc/samtools.html) |
        | Original Publication(s) | [The Sequence Alignment/Map format and SAMtools](https://doi.org/10.1093/bioinformatics/btp352)<br>[Twelve Years of SAMtools and BCFtools](https://doi.org/10.1093/gigascience/giab008) |