??? task "`extract mutation context`: Contextualizing Mpox Mutation ==(only for mpox)=="
    This task quantifies the number of G/A, C/T, and dinucleotide conversions for mpox samples. These mutations have been shown to be a characteristic of APOBEC3-type editing, which indicate adaptation of the virus to circulation among humans.

    The results are formatted into a JSON that is used to contextualize mutation metadata on the Augur phylogeny.

    When visualizing the output `auspice_input_json` file, there will be 2 new choices in the drop-down menu for "Color By":

    - G→A or C→T fraction
    - NGA/TCN context of G→A or C→T mutations.

    An example Mpox tree with these "Color By" options can be viewed [here](https://nextstrain.org/mpox/clade-IIb?c=GA_CT_fraction)

    !!! techdetails "Extract Mutation Context Technical Details"        
        |  | Links |
        | --- | --- |
        | Task | [task_augur_mutation_context.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/augur/task_augur_mutation_context.wdl) |
