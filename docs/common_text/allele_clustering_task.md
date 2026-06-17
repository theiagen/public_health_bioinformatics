---
title: Task Fragment `allele_clustering`
fragment: true
---

??? task "`allele_clustering`: PulseNet 2.0 Hash-Based Allele Clustering"
    The Allele Clustering module is used by PulseNet 2.0 to generate NWK trees for visualization, using the results from the `allele_clustering` task.

    To run this task, a tree building algorithm and distance algorithm must be specified. See available options below.

    !!! dna "`tree_building_algorithm` options"
        The options for the `tree_building_algorithm` input are as follows:

        - **"upgma"** - Iteratively merges the two closest clusters, updating inter-cluster distances with the size-weighted average of the merged cluster's distances, producing a clock-like tree. It is fast and simple but assumes a constant evolutionary rate, so it can distort topology when that assumption is violated.
        - **"single_linkage"** - Merges the two clusters whose _closest_ cross-cluster sample pair distance is the smallest across all cluster pairs, using the _minimum_ of the two merged clusters' distances. This can cause "chaining," where a sample slightly closer to one cluster than another pulls entire groups together, producing long, straggly trees
        - **"absolute_linkage"** - Merges the two clusters whose _farthest_ cross-cluster sample pair distance is the smallest across all other cluster pairs, using the _maximum_ of the two merged clusters' distances. It tends to produce compact, evenly sized clusters, but can over-split groups that are internally diverse.
        - "**neighbor_joining"** - A distance-based method that corrects raw distances by the average divergence of each sample before selecting the pair to join, producing an unrooted tree with branch lengths. It handles rate variation well, but it requires all pairwise distances up front and can produce negative branch lengths in noisy data.
        - **"minimum_spanning"** - An implementation intended to reproduce the behavior of the BioNumerics algorithm, this method builds a minimum spanning network by greedily connecting each unjoined sample to its nearest already-connected neighbor with a single edge. It produces a network rather than a bifurcating tree, making it ideal for visualizing microevolutionary relationships, but the resulting graph is not a phylogenetic tree and cannot be interpreted as one.

    !!! dna "`distance_algorithm` options"
        Distance metrics are calculated by specifying one of the following options:

        - **"absolute_allele_differences"** - Counts the raw number of loci where two samples have different alleles, ignoring any loci where either sample has a missing value. This is the simplest and most interpretable metric, but it is sensitive to the number of loci typed; samples with more loci called will naturally accumulate more differences.
        - **"normalized_allele_differences"** - Divides the count of differing alleles (including missing data) by the number of loci where _both_ samples have a call, expressing distance as a percentage (0-100). This corrects for variable typing completeness across samples, making comparisons fairer, but the percentage scale can obscure the absolute magnitude of differences.

    !!! techdetails "Allele Clustering Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_allele_clustering.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/task_allele_clustering.wdl) |
        | Software Source Code | [PulseNet 2.0 Trees](https://github.com/ncezid-biome/pulsenet2.0-trees/tree/main) |
        | Software Documentation |  [PulseNet 2.0 Trees on GitHub](https://github.com/ncezid-biome/pulsenet2.0-trees/tree/main#readme) |
