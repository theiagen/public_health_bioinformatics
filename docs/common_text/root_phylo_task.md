??? task "`root_phylo`"

    Root_Phylo returns a rooted phylogeny from inputted outgroup(s) or by finding the midpoint root. Outgroups must be tip names (case-sensitive) that exist within the tree, and multiple outgroups must be comma-delimited. Up to two outgroup tips can be supplied, and the most-recent common ancestor (MRCA) of the these outgroups will be used as the rooting branch. It is important to note that rooting on the MRCA of two outgroups is relative to the topology of the tree prior to rooting - if one of the samples is at that base of the phylogeny prior to rooting, then a random tip will be selected to allow for rooting upon the MRCA of the two inputted outgroups.

    !!! techdetails "Root_Phylo Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_root_phylo.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/utilities/task_root_phylo.wdl) |
        | Software Source Code | <https://github.com/theiagen/theiaphylo> |
        | Software Documentation | [TheiaPhylo](https://github.com/theiagen/theiaphylo/blob/main/README.md) |