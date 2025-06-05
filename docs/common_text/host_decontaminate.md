??? task "`host_decontaminate`: Host read decontamination"

    Host genetic data is frequently incidentally sequenced alongside pathogens, which can negatively affect the quality of downstream analysis. Host Decontaminate attempts to remove host reads by aligning to a reference host genome acquired on-the-fly. The reference host genome can be acquired via [NCBI Taxonomy-compatible](https://www.ncbi.nlm.nih.gov/taxonomy) taxon input or assembly accession. Host Decontaminate maps inputted reads to the host genome using `minimap2`, reports mapping statistics to this host genome, and outputs the unaligned dehosted reads. 

    !!! techdetails "Digger-Denovo Technical Details"
        |  | Links |
        | --- | --- |
        | SubWorkflow File | [wf_host_decontaminate.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/standalone/wf_host_decontaminate.wdl) |