??? task "`reorder_matrix`: Matrix Reorder and Phylogenetic Tree Rooting"

    Reorder Matrix will reorder a TSV matrix file's entries in the same order as the tips of an inputted phylogenetic tree. The phylogenetic tree can be rooted at the midpoint via the `midpoint_root_tree` Boolean, or via an outgroup tip via the `outgroup_root` String input. The provided outgroup tip must be an exact match for a sequence input into the phylogenetic tree. Please note that phylogenetic tree tip names are derived from the alignment FASTA file headers in PHB, and may deviate from `samplename` inputs.

    !!! techdetails "Reorder Matrix Technical Details"        
        |  | Links |
        | --- | --- |
        | Task | [task_reorder_matrix.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/utilities/task_reorder_matrix.wdl) |