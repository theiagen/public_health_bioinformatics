
??? task "SNP-sites (optional)"
    ##### SNP-sites (optional)

    !!! tip "Turn on SNP-Sites with `core_genome`"
        SNP-sites runs when the `core_genome` option is set to true.

    SNP-sites is used to filter out invariant sites in the whole-genome alignment, thereby creating a core genome alignment for phylogenetic inference. The output is a fasta file containing the core genome of each sample only. If Gubbins has been used, this output fasta will not contain any sites that are predicted to have arisen via recombination.

    !!! techdetails "SNP-sites technical details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_snp_sites.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/utilities/task_snp_sites.wdl) |
        | Software Source Code | [SNP-sites on GitHub](https://github.com/sanger-pathogens/snp-sites) |
        | Software Documentation | [SNP-sites on GitHub](https://github.com/sanger-pathogens/snp-sites) |
        | Original Publication(s) | [SNP-sites: rapid efficient extraction of SNPs from multi-FASTA alignments](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000056) |
