??? task "`augur traits`: Ancestral Trait Reconstruction"

    Augur Traits will reconstruct the ancestral traits of provided metadata. By default, only the "pango_lineage" and "clade_membership" columns are included, though the `augur_traits_columns` String inputted can be populated with a comma-delimited string to determine what trait metadata to use.

    !!! techdetails "Augur Traits Technical Details"        
        |  | Links |
        | --- | --- |
        | Task | [task_augur_traits.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/augur/task_augur_traits.wdl) |
        | Software Source Code | [Augur](https://github.com/nextstrain/augur) |
        | Software Documentation | [Augur traits](https://docs.nextstrain.org/projects/augur/en/stable/usage/cli/traits.html) |
        | Original Publication(s) | [Nextstrain: real-time tracking of pathogen evolution](https://doi.org/10.1093/bioinformatics/bty407) |
