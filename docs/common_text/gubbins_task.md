
??? task "Gubbins (optional)"
    ##### Gubbins (optional)

    !!! tip "Turn on Gubbins with `use_gubbins`"
        Gubbins runs when the `use_gubbins` option is set to `true` (default=true).

<!-- if: snippy_streamline -->
    **_Most optional inputs are hidden in Snippy_Streamline for simplification of the workflow. If you would like to use Gubbins with additional options, please use the `Snippy_Tree` workflow._**

    !!! tip ""
        In Snippy Streamline, the nucleotide substitution model used by gubbins will always be GTR+GAMMA.
<!-- endif -->

    **G**enealogies **U**nbiased **B**y recom**B**inations **I**n **N**ucleotide **S**equences (Gubbins) identifies and masks genomic regions that are predicted to have arisen via recombination. It works by iteratively identifying loci containing elevated densities of SNPs and constructing phylogenies based on the putative single nucleotide variants outside these regions (for more details, see [here](https://github.com/nickjcroucher/gubbins/blob/v3.3/docs/gubbins_manual.md#description-of-the-algorithm)). By default, these phylogenies are constructed using RaxML and a GTR-GAMMA nucleotide substitution model, which will be the most suitable model for most bacterial phylogenetics, though this can be modified with the `tree_builder` and `nuc_subst_model` inputs.

    Gubbins is the industry standard for masking recombination from bacterial genomes when building phylogenies, but limitations to recombination removal exist. Gubbins cannot distinguish recombination from high densities of SNPs that may result from assembly or alignment errors, mutational hotspots, or regions of the genome with relaxed selection. The tool is also intended only to find recombinant regions that are short relative to the length of the genome, so large regions of recombination may not be masked. These factors should be considered when interpreting resulting phylogenetic trees, but overwhelmingly Gubbins improves our ability to understand ancestral relationships between bacterial genomes.

<!-- if: snippy_tree -->
    There are few optional inputs for Gubbins that can be modified by the user:

    - `iterations`: Gubbins works by iteratively identifying loci containing elevated densities of SNPs, while constructing phylogenies based on the putative single nucleotide variants outside these regions. It may take many iterations for Gubbins to converge on an alignment that it considers free of recombination, especially for phylogenies that contain large numbers of genomes. By default, Gubbins is limited to 5 iterations though this may be increased by the user with the `iterations`optional input (incurring increased computing time and cost, and possibly requiring increased memory allocation).
    - `nuc_subst_model`, `tree_builder` and `tree_args`:  When Gubbins constructs phylogenies, it can use a number of phylogenetic inference tools, each with [different nucleotide substitution models](https://github.com/nickjcroucher/gubbins/blob/master/docs/gubbins_manual.md#nucleotide-substitution-model-options) and [tree-building models](https://github.com/nickjcroucher/gubbins/blob/master/docs/gubbins_manual.md#tree-building-options). By default, the `Snippy_Tree` workflow uses a GTRGAMMA substitution model and RaxML for tree building (typically suitable for bacterial genomes), but these can be modified by the user depending on the genome sequences being used with the `nuc_subst_model` and `tree_builder` optional inputs, respectively. The nucleotide substitution models that are available depend on the tree building algorithm being used (see [here](https://github.com/nickjcroucher/gubbins/blob/v3.3/docs/gubbins_manual.md#nucleotide-substitution-model-options)). Additional options for generating the phylogenetic trees in Gubbins can be specified with the `tree_args` optional input, providing an input string that is consistent with the option formats of the Gubbins command.
    - `filter_percent`: By default, Gubbins removes genomes from the multiple sequence alignment if  more than 25 % of the genome is represented by gaps. The percentage of gaps can be modified by the user using the `filter_percent` optional input.
<!-- endif -->

    !!! techdetails "Gubbins Technical Details"
                
        |  | Links |
        | --- | --- |
        | Task | [task_gubbins.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/task_gubbins.wdl) |
        | Software Source Code | [Gubbins on GitHub](https://github.com/nickjcroucher/gubbins) |
        | Software Documentation | [Gubbins v3.3 manual](https://github.com/nickjcroucher/gubbins/blob/v3.3/docs/gubbins_manual.md) |
        | Original Publication(s) | [Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences using Gubbins](https://academic.oup.com/nar/article/43/3/e15/2410982) |
    