??? task "Snippy"
    ##### Snippy
<!-- if: snippy_tree -->
    Snippy is a pipeline for calling SNPs and INDELs in haploid genomes. Before running `Snippy_Tree`, you must run `Snippy_Variants`, another workflow that uses the Snippy tool to align reads against a reference genome for individual samples. In `Snippy_Tree`, the snippy tool is used again to generate a whole-genome multiple sequence alignment (fasta file) of reads from all the samples we'd like in our tree. 
<!-- endif -->
<!-- if: snippy_streamline -->
    Snippy is used to generate a whole-genome multiple sequence alignment (fasta file) of reads from all the samples we'd like in our tree. 
<!-- endif -->

    When generating the multiple sequence alignment, a bed file can be provided by users to mask certain areas of the genome in the alignment. This is particularly relevant for masking known repetitive regions in _Mycobacterium tuberculosis_  genomes, or masking known regions containing phage sequences.

    !!! info "Why do I see `snippy_core` in Terra?"
        In Terra, this task is named "snippy_core" after the name of the command in the original Snippy tool. Despite the name, this command is NOT being used to make a core genome, but instead a multiple sequence alignment of the whole genome (without any sections masked using a bed file).
        
    !!! techdetails "Snippy Technical Details"
    
        |  | Links |
        | --- | --- |
        | Task | [task_snippy_core.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/utilities/task_snippy_core.wdl) |
        | Software Source Code | [Snippy on GitHub](https://github.com/tseemann/snippy) |
        | Software Documentation | [Snippy on GitHub](https://github.com/tseemann/snippy) |
