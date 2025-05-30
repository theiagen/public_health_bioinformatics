??? task "`KmerFinder`: Taxon Assignment (optional)"

    The `KmerFinder` method predicts prokaryotic species based on the number of overlapping (co-occurring) *k*-mers, i.e., 16-mers, between the query genome and genomes in a reference database. These *k*-mers are selected with the prefix of ATGAC as to focus on coding regions of genomes. A prediction is made by identifying which species in the training data has the highest number of 16-mers in common with the query. This match is made regardless of position. Ties will result in an alphabetical sorting of tied species and the return of the first species. This is a simpler approach than other *k*-mer based tools such as GAMBIT which utilizes a 11-mer approach while implementing thresholds to its classification algorithm allowing it to re assess at a higher level of classification if need be. 

    !!! techdetails "KmerFinder Technical Details"        
        
        |  | Links |
        | --- | --- |
        | Task | [task_kmerfinder.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/contamination/task_kmerfinder.wdl) |
        | Software Source Code | [KmerFinder BitBucket](https://bitbucket.org/genomicepidemiology/kmerfinder) |
        | Software Documentation | [CGE KmerFinder Documentation](https://cge.food.dtu.dk/services/KmerFinder/instructions.php) |
        | Original Publication(s) | [**Benchmarking of Methods for Genomic Taxonomy**](https://journals.asm.org/doi/full/10.1128/jcm.02981-13?rfr_dat=cr_pub++0pubmed&url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org) |
