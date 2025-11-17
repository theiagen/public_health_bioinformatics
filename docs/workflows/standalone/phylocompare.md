# PhyloCompare

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Name": "[**PhyloCompare**](../workflows/standalone/phylocompare.md)"}, columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level", "Dockstore"]) }}

## PhyloCompare_PHB

PhyloCompare will generate a cophylogeny plot that visualizes the differences in two trees' tip arrangements. PhyloCompare can also quantitatively compare two phylogenies by calculating the distance between two trees as a measure of the difference in their topologies (tip and branch arrangement). Validation is triggered by setting the `validate` boolean to "true".

It is recommended to root a phylogeny and PhyloCompare can root upon an outgroup tip or the midpoint.

??? dna "Tree rooting"
    If no rooting options are supplied PhyloCompare will determine if the trees are rooted or unrooted. 

    `outgroup` and `midpoint` are incompatible options and the `outgroups` input will take precedence.

??? warning "`phylovalidate_flag` errors"
    The `phylovalidate_flag` flags information that may confound distance calculation; e.g. "polytomy" can confound tree comparison if there are non-0 length branches descending from a polytomy, which may lead to erroneous distances if tips are reported in different order. In other words, phylogenies with the same topology may be reported with a non-0 distance if the tips within a polytomy are rearranged within the tree file.

    If flags are accompanied by a ">0" `phylocompare_distance`, then this indicates no distance was calculated; e.g. the "edge_count_mismatch" flag is raised when the number of edges differs between trees and a distance could not be calculated.  

### Inputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "PhyloCompare"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Workflow Tasks

{{ include_md("common_text/root_phylo_task.md", indent=0) }}

{{ include_md("common_text/cophylogeny_task.md", indent=0) }}

{{ include_md("common_text/phylovalidate_task.md", indent=0) }}

### Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "PhyloCompare"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///

## References

> - Lin, Y., Rajan, V., Moret, B. M. E. (2012). A metric for phylogenetic trees based on matching. IEEE/ACM Transactions on Computational Biology and Bioinformatics, 9(4), 1014-22, <https://doi.org/10.1109/tcbb.2011.157>

> - Moon, J. & Eulenstein, O. (2018). Cluster Matching Distance for Rooted Phylogenetic Trees. Lecture Notes in Computer Science, 10847, <https://doi.org/10.1007/978-3-319-94968-0_31>