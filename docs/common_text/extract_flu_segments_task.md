??? task "`extract_flu_segments`"

    This task takes in a full genome assembly in multifasta format containing individual Influenza segments and extracts each segment into its own fasta file. It also generates a concatenated fasta file combining all segments into one sequence.

    Labeling individual segments is based on the input flu type (Type A or Type B) and the length of each segment sequence.

    !!! techdetails "Extract Flu Segments Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_extract_flu_segments.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/influenza/task_extract_flu_segments.wdl) |