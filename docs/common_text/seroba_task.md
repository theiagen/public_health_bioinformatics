??? task "`SeroBA`: Serotyping ==_for Illumina_PE only_=="
    SeroBA is a k-mer based method for serotyping using the capsular polysaccharide biosynthesis (_cps_) locus of _Streptococcus pneumoniae_ from paired-end sequencing data. This locus encodes the serotype and is a major virulence factor for the species. Identifying circulating serotypes is important to determine the epidemiological trends and vaccine impact.

    By adapting a database from [PneumoCaT (Pneumococcal Capsular Typing)](https://github.com/ukhsa-collaboration/PneumoCaT), SeroBA uses [KMC](https://github.com/refresh-bio/KMC) to generate a k-mer database and then uses a capsular type variant database and an [ARIBA](https://github.com/sanger-pathogens/ariba)-compatible database that clusters all serotypes by serogroups. A k-mer analysis is performed and the serotype with the highest normalized sequence converage is selected. ARIBA then is used to build an assembly to confirm the selected serotype from the read data and aligns the _cps_ sequence against a reference to identify variants.

    !!! techdetails "SeroBA Technical Details"        
        |  | Links |
        | --- | --- |
        | Task | [task_seroba.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/streptococcus/task_seroba.wdl) |
        | Software Source Code | [SeroBA on GitHub](https://github.com/sanger-pathogens/seroba) |
        | Software Documentation | [SeroBA on GitHub](https://github.com/sanger-pathogens/seroba) |
        | Original Publication(s) | [SeroBA: rapid high-throughput serotyping of Streptococcus pneumoniae from whole genome sequence data](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000186) |
