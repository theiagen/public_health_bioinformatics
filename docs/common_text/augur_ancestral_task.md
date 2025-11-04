??? task "`augur ancestral`: Ancestral Nucleotide Sequence Reconstruction"
    `augur ancestral` infers ancestral nucleotide sequences based on phylogenetic relatedness using maximum-likelihood via TreeTime. A "joint" maximum likelihood model is used by default, though "marginal" input is permitted. 

    !!! techdetails "`augur ancestral` Technical Details"        
        |  | Links |
        | --- | --- |
        | Task | [task_augur_ancestral.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/augur/task_augur_ancestral.wdl) |
        | Software Source Code | [Augur on GitHub](https://github.com/nextstrain/augur)<br>[TreeTime on GitHub](https://github.com/neherlab/treetime) |
        | Software Documentation | [augur ancestral on Nexstrain](https://docs.nextstrain.org/projects/augur/en/stable/usage/cli/ancestral.html)<br>[TreeTime ReadTheDocs](https://treetime.readthedocs.io/en/latest/) |
        | Original Publication(s) | [Nextstrain: real-time tracking of pathogen evolution](https://doi.org/10.1093/bioinformatics/bty407)<br>[TreeTime: Maximum-likelihood phylodynamic analysis](https://academic.oup.com/ve/article/4/1/vex042/4794731?login=false) |
