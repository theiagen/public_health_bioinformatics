??? task "`nextclade`"

    Nextclade is an open-source project used to analyze viral genomes, particularly for clade assignment and mutation calling. Simply, Nextclade works by aligning viral genomes to a reference genome, calling variants between the two sequences, and then assigning clades based on those identified mutations. 
    
    Clade assignment is performed via _phylogenetic placement_. Phylogenetic placement compares the mutations of the provided sequence to the mutations of each node found in a reference tree, where the root of that tree is the reference genome. The node that is most similar to the sample is used to both assign a clade designation and calculate where the sample should be placed in the phylogenetic tree.
    
<!-- if: theiaviral -->
    Theiagen has implemented a full genome-based [Nextclade dataset](https://github.com/theiagen/rabies) for *L. rabies* with subclade classification resolution.
<!-- endif -->

    !!! techdetails "Nextclade Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_nextclade.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/task_nextclade.wdl) |
        | Software Source Code | <https://github.com/nextstrain/nextclade> |
        | Software Documentation | [Nextclade](https://docs.nextstrain.org/projects/nextclade/en/stable/) |
        | Original Publication(s) | [Nextclade: clade assignment, mutation calling and quality control for viral genomes.](https://doi.org/10.21105/joss.03773) |
