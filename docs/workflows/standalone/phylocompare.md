# PhyloCompare

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Standalone](../../workflows_overview/workflows_type.md/#standalone) | [Any Taxa](../../workflows_overview/workflows_kingdom.md/#any-taxa) | PHB vX.X.X | No | Sample |

## PhyloCompare_PHB

PhyloCompare will calculate the distance between two _newick-formatted_ phylogenies as a measure of the difference in their topologies (tip and branch arrangement). Additionally, PhyloCompare will validate if the phylogenies exceed an inputted maximum distance. The maximum distance is 0 by default, which indicates that the phylogenies have the same topology. 

It is difficult to conceptualize what a non-0 distance indicates, so please see the following citations for their interpretation. For unrooted phylogenies, PhyloCompare calculates the [Lin-Rajan-Moret distance](https://pubmed.ncbi.nlm.nih.gov/22184263/), and for rooted phylogenies, PhyloCompare calculates the [matching cluster distance](https://link.springer.com/chapter/10.1007/978-3-319-94968-0_31#:~:text=Phylogenetic%20trees%20are%20fundamental%20to%20biology%20and,is%20an%20important%20problem%20in%20computational%20phylogenetics.). The Robinson-Foulds distance is also calculated, though it is disregarded in validation (see citations for criticism).

PhyloCompare can automatically root upon outgroup tips or the midpoint. If more than a single outgroup tip is supplied then the phylogenies will be rooted on their most recent common ancestor branch. Input multiple outgroup tips as a comma-delimited list, e.g. "tip1,tip2". 

!!! warning
    If no rooting options are supplied and `unrooted` is not set to `true`, PhyloCompare will assume the supplied trees are rooted. 

    `unrooted`, `root_tips`, and `midpoint` are incompatible options and selecting multiple in the same run may lead to undesired behavior.

### Inputs

<div class="searchable-table" markdown="1">

Please note that all string inputs **must** be enclosed in quotation marks; for example, "tip1,tip2" or "tip1".

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| phylocompare | **tree1_path** | String | Path to a newick-formatted phylogenetic tree in an accessible bucket |  | Required |
| phylocompare | **tree2_path** | String | Path to a newick-formatted phylogenetic tree in an accessible bucket |  | Required |
| phylocompare | **max_distance** | Float | Maximum tolerable distance in validation | 0.0 | Optional |
| phylocompare | **midpoint** | Boolean | Root phylogenies at their midpoint | false | Optional |
| phylocompare | **root_tips** | String | Comma-delimited list of outgroup tip(s) to root upon. Multiple outgroup tips will root on the branch descended from their most recent common ancestor | | Optional |
| phylocompare | **unrooted** | Boolean | Interpret the phylogenies as unrooted | false | Optional |
| phylovalidate | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional |
| phylovalidate | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| phylovalidate | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/theiaphylo:0.1.0 | Optional |
| phylovalidate | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 4 | Optional |
| version_capture | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional |

</div>

### Outputs

<div class="searchable-table" markdown="1">

| **Variable** | **Type** | **Description** |
|---|---|---|
| phb_version | String | PHB version |
| phylo_distance | String | Distance between the phylogenies, or "None" if distance was unable to be calculated |
| phylocompare_report | File | Text file of the calculated distances |
| phylocompare_version | String | Version of PhyloCompare python script |
| validation | String | "PASS" if distance < `max_distance` and "FAIL" if distance > `max_distance` or could not be calculated |

</div>

## References

> - Lin, Y., Rajan, V., Moret, B. M. E. (2012). A metric for phylogenetic trees based on matching. IEEE/ACM Transactions on Computational Biology and Bioinformatics, 9(4), 1014-22, <https://doi.org/10.1109/tcbb.2011.157>

> - Moon, J. & Eulenstein, O. (2018). Cluster Matching Distance for Rooted Phylogenetic Trees. Lecture Notes in Computer Science, 10847, <https://doi.org/10.1007/978-3-319-94968-0_31>

> - Cogent3 Python Library <https://github.com/cogent3/cogent3>
