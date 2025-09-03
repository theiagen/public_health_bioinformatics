??? task "`PlasmidFinder`: Plasmid Identification (optional)"

    Activate this task by setting `call_plasmidfinder` to `true`.

    `PlasmidFinder`, from the Center for Genomic Epidemiology, detects plasmids in total or partially sequenced genomes and identifies the closest plasmid type in the database for typing purposes.

    ???+ toggle "What are plasmids?"
        
        Plasmids are double-stranded circular or linear DNA molecules that are capable of replication independently of the chromosome and may be transferred between different species and clones. Many plasmids contain resistance or virulence genes, though some do not clearly confer an advantage to their host bacterium.
        
    !!! techdetails "PlasmidFinder Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_plasmidfinder.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/plasmid_detection/task_plasmidfinder.wdl) |
        | Software Source Code | [PlasmidFinder Tool on BitBucket](https://bitbucket.org/genomicepidemiology/plasmidfinder/src/master/)<br>[PlasmidFinder Database on BitBucket](https://bitbucket.org/genomicepidemiology/plasmidfinder_db/src/master/) |
        | Software Documentation | [PlasmidFinder on BitBucket](https://bitbucket.org/genomicepidemiology/plasmidfinder/src/master/) |
        | Original Publication(s) | [In Silico Detection and Typing of Plasmids using PlasmidFinder and Plasmid Multilocus Sequence Typing](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4068535/)<br> |
