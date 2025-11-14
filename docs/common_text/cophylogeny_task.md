??? task "`cophylogeny`"

    The Cophylogeny task will generate cophylogeny plots of two inputted phylogenies. A cophylogeny plot draws lines between two trees' tips as a method for visualizing their topological (tip/branch arrangment) differences. 

    A cophylogeny plot is generated with branch lengths (`cophylogeny_plot_with_branch_lengths`) and a cophylogeny plot without branch lengths (`cophylogeny_plot`). The plot without branch lengths is better for depicting branching order differences, though it is important to note that the branch lengths within this plot are arbitrary and do not convey evolutionary distance. Users will most likely need to visualize the phylogenies independently to interpret evolutionary distance because it is difficult to automatically graph two phylogenies with scaled and viewable branch lengths.

    !!! techdetails "Cophylogeny Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_cophylogeny.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/utilities/task_cophylogeny.wdl) |
        | Software Source Code | <https://github.com/theiagen/theiaphylo> |
        | Software Documentation | [TheiaPhylo](https://github.com/theiagen/theiaphylo/blob/main/README.md) |