??? task "`PopPUNK`: Global Pneumococcal Sequence Cluster Typing"
    Global Pneumococcal Sequence Clusters (GPSC) define and name pneumococcal strains. Each GPSC is an [international definition of a pneumococcal lineage](https://doi.org/10.1016/j.ebiom.2019.04.021) and can capture all variations across the entire genome, leading to better vaccine development due to increased disease surveillance. GPSC designation is undertaken using the PopPUNK (Population Partitioning Using Nucleotide K-mers) software and the GPSC database.

    PopPUNK works by using variable-length k-mers to distinguish between sample divergences in shared genomic content. Clusters are assigned based on the resulting pairwise distance distributions. 

    !!! tip "Interpreting GPSC results"
        - In the `*_external_clusters.csv` novel clusters are assigned NA. For isolates that are assigned a novel cluster and pass QC, you can email [globalpneumoseq@gmail.com](mailto:globalpneumoseq@gmail.com) to have these novel clusters added to the database.
        - Unsampled diversity in the pneumococcal population may represent missing variation that links two GPS clusters. When this is discovered, GPSCs are merged and the merge history is indicated. For example, if GPSC23 and GPSC362 merged, the GPSC would be reported as GPSC23, with a merge history of GPSC23;362.
    
    !!! techdetails "PopPUNK Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_poppunk_streppneumo.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/streptococcus/task_poppunk_streppneumo.wdl) |
        | Software Source Code | [PopPUNK on GitHub](https://github.com/bacpop/PopPUNK)<br>[Global Pneumococcal Sequencing Project](https://www.pneumogen.net/gps/#/)<br>[GPSC Database](https://www.pneumogen.net/gps/#/gpsc#lineages) |
        | Software Documentation | [PopPUNK Documentation](https://poppunk-docs.bacpop.org/)<br>[GPS Training Documentation](https://www.pneumogen.net/gps/#/training) |
        | Original Publication(s) | _PopPUNK tool_: [Fast and flexible bacterial genomic epidemiology with PopPUNK](https://genome.cshlp.org/content/29/2/304)<br>_GPSC database_: [International genomic definition of pneumococcal lineages, to contextualise disease, antibiotic resistance and vaccine impact](https://doi.org/10.1016/j.ebiom.2019.04.021) |
    