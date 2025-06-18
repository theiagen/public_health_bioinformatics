??? task "`host_decontaminate`: Host read decontamination"

    Host genetic data is frequently incidentally sequenced alongside pathogens, which can negatively affect the quality of downstream analysis. Host Decontaminate attempts to remove host reads by aligning to a reference host genome acquired on-the-fly. The reference host genome can be acquired via [NCBI Taxonomy-compatible](https://www.ncbi.nlm.nih.gov/taxonomy) taxon input or assembly accession. Host Decontaminate maps inputted reads to the host genome using `minimap2`, reports mapping statistics to this host genome, and outputs the unaligned dehosted reads. 

    A common alternative approach to decontamination is to use a fast, usually k-mer-based, read classification software such as `kraken2`/`Metabuli`. These software have the advantage of screening large databases of many organisms much faster than full read alignment, though the quality of decontamination is confounded by high rates of false negative/false positive alignments compared to direct alignment. When the host is known, direct alignment to a reference host genome via `minimap2` is likely to yield higher quality mapping, which is the underlying justification for Host Decontaminate's approach.

    The detailed steps and tasks are as follows:

    ??? toggle "Taxonomic Identification"

        {{ include_md("common_text/ncbi_identify_task.md", condition="host_decontaminate", indent=8) }}

    ??? toggle "Download Accession"

        {{ include_md("common_text/ncbi_identify_task.md", condition="host_decontaminate", indent=8) }}

    ??? Toggle "Map Reads to Host"

        {{ include_md("common_text/minimap2_task.md", condition="host_decontaminate", indent=8) }}

    



    !!! techdetails "Host Decontaminate Technical Details"
        |  | Links |
        | --- | --- |
        | SubWorkflow File | [wf_host_decontaminate.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/standalone/wf_host_decontaminate.wdl) |