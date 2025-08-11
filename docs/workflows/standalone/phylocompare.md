# PhyloCompare

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Name": "[**PhyloCompare**](../workflows/standalone/phylocompare.md)"}, columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level", "Dockstore"]) }}

## PhyloCompare_PHB

PhyloCompare will calculate the distance between two _newick-formatted_ phylogenies as a measure of the difference in their topologies (tip and branch arrangement). A distance of 0 indicates the phylogenies have the same topology. PhyloCompare will validate if the phylogenies exceed an inputted maximum distance. The maximum distance is 0 by default and the phylogenies must have the same tips.

It is difficult to conceptualize what a non-0 distance indicates, so please see the following citations for their interpretation. For unrooted phylogenies, PhyloCompare calculates the [Lin-Rajan-Moret distance](https://pubmed.ncbi.nlm.nih.gov/22184263/), and for rooted phylogenies, PhyloCompare calculates the [matching cluster distance](https://link.springer.com/chapter/10.1007/978-3-319-94968-0_31#:~:text=Phylogenetic%20trees%20are%20fundamental%20to%20biology%20and,is%20an%20important%20problem%20in%20computational%20phylogenetics.). The Robinson-Foulds distance is also calculated, though it is disregarded in validation (see citations for criticism).

PhyloCompare can automatically root upon outgroup tips or the midpoint. If more than a single outgroup tip is supplied then the phylogenies will be rooted on their most recent common ancestor branch. Input multiple outgroup tips as a comma-delimited list, e.g. "tip1,tip2". 

??? dna "Tree rooting"
    If no rooting options are supplied PhyloCompare will determine if the trees are rooted or unrooted. 

    `outgroups` and `midpoint` are incompatible options and the `outgroups` input will take precedence.

??? warning "`phylocompare_flag` errors"
    The `phylocompare_flag` flags information that may confound distance calculation; e.g. "polytomy" can confound tree comparison if there are non-0 length branches descending from a polytomy, which may lead to erroneous distances if tips are reported in different order. In other words, phylogenies with the same topology may be reported with a non-0 distance if the tips within a polytomy are rearranged within the tree file.

    If flags are accompanied by a ">0" `phylocompare_distance`, then this indicates no distance was calculated; e.g. the "edge_count_mismatch" flag is raised when the number of edges differs between trees and a distance could not be calculated.  

### Inputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "PhyloCompare"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "PhyloCompare"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///

### All Tasks

{{ include_md("common_text/root_phylo_task.md", indent=0) }}

{{ include_md("common_text/phylocompare_task.md", indent=0) }}

## References

> - Lin, Y., Rajan, V., Moret, B. M. E. (2012). A metric for phylogenetic trees based on matching. IEEE/ACM Transactions on Computational Biology and Bioinformatics, 9(4), 1014-22, <https://doi.org/10.1109/tcbb.2011.157>

> - Moon, J. & Eulenstein, O. (2018). Cluster Matching Distance for Rooted Phylogenetic Trees. Lecture Notes in Computer Science, 10847, <https://doi.org/10.1007/978-3-319-94968-0_31>

> - Cogent3 Python Library <https://github.com/cogent3/cogent3>
