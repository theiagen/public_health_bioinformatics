---
title: Task Fragment `template`
fragment: true
---
??? task "`PIRATE`: Gene Family Clustering and Alignment"
    PIRATE uses GFF3 files to classify genes into orthologous gene families in bacterial pangenomes using sequence identity thresholds. By default, genes are aligned and a core/pangenome alignment is produced. This option can be turned off by setting `align` to "false". The pangenome is constructed using the CDS (modifiable with the `features` parameter) and are not translated to amino acid sequences (modifiable with the `nucl` parameter). 

    PIRATE generates (by default) a pangenome summary (the number and frequency of enes in the pangenome); tabular summaries of the gene families and the unique alleles belonging to each family; a binary tree using presence-absence data; and core and pangenome alignment files and their respective annotation files. 
  
    !!! techdetails "PIRATE Technical Details"        
        |  | Links |
        | --- | --- |
        | Task | [task_pirate.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/task_pirate.wdl) |
        | Software Source Code | [PIRATE on GitHub](https://github.com/SionBayliss/PIRATE) |
        | Software Documentation |  [PIRATE on GitHub](https://github.com/SionBayliss/PIRATE) |
        | Original Publication(s) | [PIRATE: A fast and scalable pangenomics toolbox for clustering diverged orthologues in bacteria>](https://doi.org/10.1093/gigascience/giz119) |
