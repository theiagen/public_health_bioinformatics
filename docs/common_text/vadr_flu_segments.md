??? task "`vadr_flu_segments`"

    This task processes a full or partial influenza genome assembly in multifasta format, along with the output `.tar.gz` file from a VADR run. It extracts each segment into its own fasta file and also generates a concatenated fasta containing all segments combined into a single sequence. Segment names are assigned based on the specified flu type (A or B) and the segment classification found in the [VADR .sqc file](https://github.com/ncbi/vadr/blob/master/documentation/formats.md#sqc).

    Note: Results may be unreliable if segment lengths deviate from those expected for Influenza A or B. For best results, the input assembly should contain all 8 segments as separate contigs.

    !!! techdetails "VADR Flu Segments Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_vadr_flu_segments.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/influenza/task_vadr_flu_segments.wdl) |