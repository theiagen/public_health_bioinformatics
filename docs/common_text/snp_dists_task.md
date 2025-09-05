??? task "SNP-dists"
<!-- if: snippy -->  
    ##### SNP-dists
<!-- endif -->
    `SNP-dists` computes pairwise SNP distances between genomes. It takes the same alignment of genomes used to generate your phylogenetic tree and produces a matrix of pairwise SNP distances between sequences. This means that if you generated pairwise core-genome phylogeny, the output will consist of pairwise core-genome SNP (cgSNP) distances. Otherwise, these will be whole-genome SNP distances. Regardless of whether core-genome or whole-genome SNPs, this SNP distance matrix will exclude all SNPs in masked regions (i.e. masked with a bed file or gubbins). 

    The SNP-distance output can be visualized using software such as [Phandango](http://jameshadfield.github.io/phandango/#/main) to explore the relationships between the genomic sequences. The task can optionally add a Phandango coloring tag (:c1) to the column names in the output matrix to ensure that all columns are colored with the same color scheme throughout by setting `phandango_coloring` to `true`.

    !!! techdetails "SNP-dists Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_snp_dists.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/utilities/task_snp_dists.wdl) |
        | Software Source Code | [SNP-dists on GitHub](https://github.com/tseemann/snp-dists) |
        | Software Documentation | [SNP-dists on GitHub](https://github.com/tseemann/snp-dists) |
