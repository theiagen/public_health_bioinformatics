
??? task "IQTree2"
    ##### IQTree2

    IQTree2 is used to build the final phylogeny. It uses the alignment generated in the previous steps of the workflow. The contents of this alignment will depend on whether any sites were masked with recombination.

    The phylogeny is generated using the maximum-likelihood method and a specified nucleotide substitution model. By default, the Snippy_Tree workflow will run Model Finder to determine the most appropriate nucleotide substitution model for your data, but you may specify the nucleotide substitution model yourself using the `iqtree2_model` optional input (see [here](http://www.iqtree.org/doc/Substitution-Models) for available models).

    IQTree will perform assessments of the tree using the Shimodaira–Hasegawa approximate likelihood-ratio test ([SH-aLRT test](https://academic.oup.com/sysbio/article/59/3/307/1702850?login=false)), and ultrafast bootstrapping with [UFBoot2](https://academic.oup.com/mbe/article/35/2/518/4565479), a quicker but less biased alternative to standard bootstrapping. A clade should not typically be trusted if it has less than 80% support from the SH-aLRT test and less than 95% support with ultrafast bootstrapping.

    !!! tip "Nucleotide substitution model"
        When `core_genome`= `true`, the default nucleotide substitution model is set to the General Time Reverside model with Gamma distribution (GTR+G). 
        
        When the user sets `core_genome`= `false`, the default nucleotide substitution model is set to the General Time Reversible model with invariant sites and Gamma distribution (`GTR+I+G`).
                
    !!! techdetails "IQTree2 technical details"
                    
        |  | Links |
        | --- | --- |
        | Task | [task_iqtree2.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/task_iqtree2.wdl) |
        | Software Source Code | [IQ-TREE on GitHub](https://github.com/iqtree/iqtree2) |
        | Software Documentation | [IQTree documentation](http://www.iqtree.org/doc/) for the latest version (not necessarily the version used in this workflow) |
        | Original Publication(s) | [IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era](https://academic.oup.com/mbe/article/37/5/1530/5721363)<br>[New Algorithms and Methods to Estimate Maximum-Likelihood Phylogenies: Assessing the Performance of PhyML 3.0](https://academic.oup.com/sysbio/article/59/3/307/1702850?login=false)<br>[Ultrafast Approximation for Phylogenetic Bootstrap](https://academic.oup.com/mbe/article/30/5/1188/997508?login=false)<br> [UFBoot2: Improving the Ultrafast Bootstrap Approximation](https://academic.oup.com/mbe/article/35/2/518/4565479?login=false)<br>[ModelFinder: fast model selection for accurate phylogenetic estimates](https://www.nature.com/articles/nmeth.4285) |
