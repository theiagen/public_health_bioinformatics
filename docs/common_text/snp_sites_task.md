---
title: Task Fragment `snp_sites`
fragment: true
---
??? task "`SNP-sites`: Genome Filtering"
<!-- if: snippy -->
    ##### SNP-sites (optional)
    !!! tip "Turn on SNP-sites with `core_genome`"
        SNP-sites runs when the `core_genome` option is set to true.

    If Gubbins has been used, the output file will not contain any sites that are predicted to have arisen via recombination.
<!-- endif -->
    
    SNP-sites is used to identify variants in a multi-FASTA alignment, and returns _only_ the sites with SNPs in FASTA format. 

    For example, if your input FASTA is as follows:

    ```
    >sample1
    ACGT
    >sample2
    ACCT
    ```

    SNP-sites will identify only the SNPs, like so:

    ```
    >sample1
    G
    >sample2
    C
    ```

    By default, wildcard bases are output. To turn off this behavior and restrict the output to only ACGT bases, set the `allow_wildcard_bases` parameter to `"false"`.

    !!! techdetails "SNP-sites technical details"
        |  | Links |
        | --- | --- |
        | Task | [task_snp_sites.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/utilities/task_snp_sites.wdl) |
        | Software Source Code | [SNP-sites on GitHub](https://github.com/sanger-pathogens/snp-sites) |
        | Software Documentation | [SNP-sites on GitHub](https://github.com/sanger-pathogens/snp-sites) |
        | Original Publication(s) | [SNP-sites: rapid efficient extraction of SNPs from multi-FASTA alignments](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000056) |
