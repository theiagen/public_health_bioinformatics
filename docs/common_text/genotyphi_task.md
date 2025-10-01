??? task "`genotyphi`: _Salmonella_ Typhi Characterization ==_for Illumina and ONT only_=="
    [`genotyphi`](https://github.com/typhoidgenomics/genotyphi) is activated upon the identification of the "Typhi" serotype by SISTR or SeqSero2 (via either the `seqsero2_predicted_serotype`, or the `sistr_predicted_serotype`). `genotyphi` divides the _Salmonella enterica_ serovar Typhi population into detailed lineages, clades, and subclades. It also detects mutations in the quinolone-resistance determining regions, acquired antimicrobial resistance genes, plasmid replicons, and subtypes of the IncHI1 plasmid which is associated with multidrug resistance.

    This task uses [Mykrobe](https://github.com/typhoidgenomics/genotyphi/blob/main/README.md#typhi-mykrobe) in order to perform k-mer based genotyping with a genotyping scheme specific to _Salmonella_ Typhi, and then parses that output using the `genotyphi` tool. This scheme divides the _Salmonella_ Typhi population into genotypes based on unique SNV markers.

    !!! techdetails "genotyphi Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_genotyphi.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/salmonella/task_genotyphi.wdl) |
        | Software Source Code | [genotyphi on GitHub](https://github.com/typhoidgenomics/genotyphi) |
        | Software Documentation | [genotyphi on GitHub](https://github.com/typhoidgenomics/genotyphi/master/README.md) |
        | Orginal publication(s) | [An extended genotyping framework for Salmonella enterica serovar Typhi, the cause of human typhoid](https://www.nature.com/articles/ncomms12827/)<br>[Five Years of GenoTyphi: Updates to the Global Salmonella Typhi Genotyping Framework](https://academic.oup.com/jid/article/224/Supplement_7/S775/6358992?login=false)<br>[Typhi Mykrobe: fast and accurate lineage identification and antimicrobial resistance genotyping directly from sequence reads for the typhoid fever agent _Salmonella_ Typhi](https://www.biorxiv.org/content/10.1101/2024.09.30.613582v1) |
