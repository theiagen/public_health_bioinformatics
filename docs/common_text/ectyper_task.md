    ??? task "`ECTyper`: Serotyping"
        
        [ECTyper](https://github.com/phac-nml/ecoli_serotyping) is a serotyping module for *E. coli*. In TheiaProk, we are using assembly files as input. ECTyper provides species identification and quality control for *E. coli* allowiung for complete reports on serotyping, Shiga toxin typing, and pathotyping. Pathotype is identified using an ECTyper internal typing database that looks at toxin and pathotype signature marker sequences. Pathotypes are used to group *E. coli* specimens based on their identified pathogenicity. 
        
        !!! techdetails "ECTyper Technical Details"
            
            |  | Links |
            | --- | --- |
            | Task | [task_ectyper.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/escherichia_shigella/task_ectyper.wdl) |
            | Software Source Code | [ECTyper on GitHub](https://github.com/phac-nml/ecoli_serotyping) |
            | Software Documentation | [ECTyper on GitHub](https://github.com/phac-nml/ecoli_serotyping) |
            | Orginal publication | [ECTyper: in silico Escherichia coli serotype and species prediction from raw and assembled whole-genome sequence data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8767331/) |