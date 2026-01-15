??? task "`VirulenceFinder`: Virulence Gene Identification"
    VirulenceFinder, from the Center for Genomic Epidemiology (GCE) in TheiaProk is only run on assembly files due to issues regarding discordant results when using read files on the web application versus the command-line. VirulenceFinder uses BLAST and KMA to query against a database of virulence genes in _E. coli_ to identify any virulence factors in the sample.

    !!! techdetails "VirulenceFinder Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_virulencefinder.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/escherichia_shigella/task_virulencefinder.wdl) |
        | Software Source Code | [VirulenceFinder on BitBucket](https://bitbucket.org/genomicepidemiology/virulencefinder/src/master/)<br>[VirulenceFinder Database on BitBucket](https://bitbucket.org/genomicepidemiology/virulencefinder_db/src/master/) |
        | Software Documentation | [VirulenceFinder on BitBucket](https://bitbucket.org/genomicepidemiology/virulencefinder/src/master/)<br>[VirulenceFinder Database on BitBucket](https://bitbucket.org/genomicepidemiology/virulencefinder_db/src/master/) |
        | Original Publication(s) | [Real-time whole-genome sequencing for routine typing, surveillance, and outbreak detection of verotoxigenic Escherichia coli](https://pubmed.ncbi.nlm.nih.gov/24574290/) |
