??? task "`phylocompare`"

    PhyloCompare will clean two phylogenies and validate if the distance between these two phylogenies' topologies is less than an inputted `max_distance` float (0 by default). Phylogenies are cleaned by converting 0 branch length nodes into polytomies, and any detected polytomies are reported as a flag. Polytomies may arbitrarily yield a non-0 distance, though if a 0 distance is reported with a polytomy then it indicates that the polytomy did not confound distance calculation.

    For unrooted phylogenies, PhyloCompare calculates the [Lin-Rajan-Moret distance](https://pubmed.ncbi.nlm.nih.gov/22184263/), and for rooted phylogenies, PhyloCompare calculates the [matching cluster distance](https://link.springer.com/chapter/10.1007/978-3-319-94968-0_31#:~:text=Phylogenetic%20trees%20are%20fundamental%20to%20biology%20and,is%20an%20important%20problem%20in%20computational%20phylogenetics.). The Robinson-Foulds distance is also calculated, though it is disregarded in validation (see citations for criticism).

    !!! techdetails "PhyloCompare Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_root_phylo.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/utilities/task_phylocompare.wdl) |
        | Software Source Code | <https://github.com/theiagen/theiaphylo> |
        | Software Documentation | [TheiaPhylo](https://github.com/theiagen/theiaphylo/blob/main/README.md) |