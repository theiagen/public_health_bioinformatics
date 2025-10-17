??? task "`augur clades`: Assigning Clades based on Sequence Signatures"

    Augur Clades assigns clades to nodes in a tree based on amino-acid or nucleotide signatures. A `clades_tsv` is required to delineate clade-defining mutations in a tab-delimitted file with the header "clade\tgene\tsite\alt", where the site is an integer site and the alternative sequence character is delineated in the "alt" column. Clade designation preferentially references amino acid sequences if they are provided.

    !!! techdetails "Augur clades Technical Details"        
        |  | Links |
        | --- | --- |
        | Task | [task_augur_clades.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/augur/task_augur_clades.wdl) |
        | Software Source Code | [Augur](https://github.com/nextstrain/augur)<br>[TreeTime](https://github.com/neherlab/treetime) |
        | Software Documentation | [Augur clades](https://docs.nextstrain.org/projects/augur/en/stable/usage/cli/clades.html) |
        | Original Publication(s) | [Nextstrain: real-time tracking of pathogen evolution](https://doi.org/10.1093/bioinformatics/bty407) |
