??? task "`extract clade mutations`: Extract Clade-Defining Signature Sequences"

    Extract Clade Mutations will create an Augur-compatible "clades.tsv" by extracting signature clade-defining sequences. A nucleotide JSON outputted by Augur Ancestral is required, and an optional amino acid JSON outputted by Augur Translate can be used to infer specific amino acid mutations. 

    Clade-defining signatures can only be extracted from monophyletic clades with unique mutation signatures. If no clade-defining mutations are reported, an error is raised. 

    !!! techdetails "Augur clades Technical Details"        
        |  | Links |
        | --- | --- |
        | Task | [task_extract_clade_mutations.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/augur/task_extract_clade_mutations.wdl) |
        | Software Source Code | [Theiagen Utilities](https://github.com/theiagen/utilities)<br>[TheiaPhylo](https://github.com/theiagen/theiaphylo) |
