---
title: Task Fragment `allele_clustering`
fragment: true
---

??? task "`allele_clustering`: PulseNet 2.0 Hash-Based Allele Clustering"
    The Allele Clustering module is used by PulseNet 2.0 to generate NWK trees for visualization, using the results from the `allele_clustering` task.

    To run this task, a tree building algorithm and distance algorithm must be specified; these options are available in the [inputs](#inputs) section of the workflow documentation.

    !!! techdetails "Allele Clustering Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_allele_clustering.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/task_allele_clustering.wdl) |
        | Software Source Code | [PulseNet 2.0 Trees](https://github.com/ncezid-biome/pulsenet2.0-trees/tree/main) |
        | Software Documentation |  [PulseNet 2.0 Trees on GitHub](https://github.com/ncezid-biome/pulsenet2.0-trees/tree/main#readme) |
