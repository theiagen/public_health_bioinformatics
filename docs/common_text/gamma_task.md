??? task "`GAMMA`: Genotyping (optional)"

    [GAMMA](https://github.com/rastanton/GAMMA/tree/main) (Gene Allele Mutation Microbial Assessment) is a protein identity based tool that identifies gene matches in microbial genomic data. GAMMA will also translate and annotate each match providing mutational and truncation information for identified matches. This is done much like AMRFinder, however, GAMMA uses BLAT as opposed to BLAST allowing for faster calls with comparable accuracy.  

    GAMMA utilizes a multifasta database of the coding sequences of genes specified by the user. This allows for GAMMA to search for AMR, hypervirulence, plasmid markers, or any prokaryotic database. The default for `task_gamma.wdl` is the GAMMA provided [ResFinder Database](gs://theiagen-public-resources-rp/reference_data/databases/gamma/default_ResFinderDB_Combined_05-06-20.fsa). GAMMA will then return gene matches with mutation and truncation information as a `.gamma` file which can be accompanied with a GFF output utilizing the `--gff` flag.

    GAMMA also allows for the usage of only nucleotide sequences rather than translated sequences. Using GAMMA-S, enabled with the boolean `run_gammas`, will find the best matches from a multifasta database without translating sequences. 

    !!! techdetails "AMRFinderPlus Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_gamma.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/awh-oh-integration-updates-dev/tasks/gene_typing/drug_resistance/task_gamma.wdl) |
        | Software Source Code | [GAMMA on GitHub](https://github.com/rastanton/GAMMA/tree/main) |
        | Software Documentation | https://github.com/rastanton/GAMMA/tree/main |
        | Original Publication(s) | [GAMMA: a tool for the rapid identification, classification and annotation of translated gene matches from sequencing data](https://academic.oup.com/bioinformatics/article/38/2/546/6355578) |