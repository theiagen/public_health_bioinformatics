
??? task "Centroid"
    ##### Centroid
    
    Centroid selects the most central genome among a list of assemblies by computing pairwise mash distances. In `Snippy_Streamline`, this centroid assembly is then used to find a closely related reference genome that can be used to generate the tree.  In order to use `Centroid`, should complete the `samplenames` input. 

    !!! techdetails "`centroid` Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_centroid.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/utilities/task_centroid.wdl) |
        | Software Source Code | <https://github.com/theiagen/centroid> |
        | Software Documentation | <https://github.com/theiagen/centroid> |
