??? task "`host_decontaminate`: Host read decontamination"

    Host genetic data is frequently incidentally sequenced alongside pathogens, which can negatively affect the quality of downstream analysis. Host Decontaminate attempts to remove host reads by aligning to a reference host genome acquired on-the-fly. The reference host genome can be acquired via [NCBI Taxonomy-compatible](https://www.ncbi.nlm.nih.gov/taxonomy) taxon input or assembly accession. Host Decontaminate maps inputted reads to the host genome using `minimap2`, reports mapping statistics to this host genome, and outputs the unaligned dehosted reads. 

    The detailed steps and tasks are as follows:

    ??? toggle "Taxonomic Identification"

{{ include_md("common_text/ncbi_identify_task.md", condition="host_decontaminate", indent=8) }}

    ??? toggle "Download Accession"

{{ include_md("common_text/ncbi_identify_task.md", condition="host_decontaminate", indent=8) }}

    ??? Toggle "Map Reads to Host"

{{ include_md("common_text/minimap2_task.md", condition="host_decontaminate", indent=8) }}

    ??? Toggle "Extract Unaligned Reads"

{{ include_md("common_text/parse_mapping_task.md", condition="sam_to_sorted_bam", indent=8) }}

{{ include_md("common_text/parse_mapping_task.md", condition="bam_to_unaligned_fastq", indent=8) }}
    
    ??? Toggle "Host Read Mapping Statistics"

{{ include_md("common_text/assembly_metrics_task.md", indent=8)}}

    !!! techdetails "Host Decontaminate Technical Details"
        |  | Links |
        | --- | --- |
        | SubWorkflow File | [wf_host_decontaminate.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_host_decontaminate.wdl) |