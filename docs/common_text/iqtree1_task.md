---
title: Task Fragment `IQ-TREE 2`
fragment: true
---
??? task "`IQ-TREE`: Phylogeny Construction"
    IQ-TREE is used to build a phylogeny using the maximum-likelihood method and a specified nucleotide substitution model.
    
    IQ-TREE will perform assessments of the tree using the Shimodaira–Hasegawa approximate likelihood-ratio test ([SH-aLRT test](https://academic.oup.com/sysbio/article/59/3/307/1702850?login=false)), and ultrafast bootstrapping with [UFBoot2](https://academic.oup.com/mbe/article/35/2/518/4565479), a quicker but less biased alternative to standard bootstrapping. A clade should not typically be trusted if it has less than 80% support from the SH-aLRT test and less than 95% support with ultrafast bootstrapping.

    !!! tip "Nucleotide substitution model"
        The default nucleotide substitution model is set to the General Time Reverside model with invariant sites and Gamma distribution (GTR+I+G). 
                
    !!! techdetails "IQ-TREE technical details"
        |  | Links |
        | --- | --- |
        | Task | [task_iqtree.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/task_iqtree.wdl) |
        | Software Source Code | [IQ-TREE on GitHub](https://github.com/iqtree/iqtree2/tree/v1.6.7) |
        | Software Documentation | [IQ-TREE on GitHub](https://github.com/iqtree/iqtree2/tree/v1.6.7) |
        | Original Publication(s) | [IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies](https://doi.org/10.1093/molbev/msu300)<br>[UFBoot2: Improving the Ultrafast Bootstrap Approximation](https://doi.org/10.1093/molbev/msx281) |
