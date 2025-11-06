??? task "`host_decontaminate`: Host Read Decontamination"

    Host genetic data is frequently incidentally sequenced alongside pathogens, which can negatively affect the quality of downstream analysis. Host Decontaminate attempts to remove host reads by aligning to a reference host genome that is directly inputted or acquired on-the-fly. The reference host genome can be inputted into the `host` input field as an assembly file (with `is_genome` set to "true"), acquired via [NCBI Taxonomy-compatible](https://www.ncbi.nlm.nih.gov/taxonomy) taxon input, or assembly accession (with `is_accession` set to "true"). Host Decontaminate maps inputted reads to the host genome using `minimap2`, reports mapping statistics to this host genome, and outputs the unaligned dehosted reads. 

    The detailed steps and tasks are as follows:

{{ include_md("common_text/ncbi_identify_task.md", indent=4, replacements={'??? task "`ncbi_identify`"' : '??? toggle "Taxonomic Identification"'}) }}

{{ include_md("common_text/ncbi_datasets_task.md", condition="theiaviral", indent=4, replacements={'??? task "NCBI Datasets"' : '??? toggle "Download Accession"'}) }}

{{ include_md("common_text/minimap2_task.md", condition="only_map_ont", indent=4, replacements={'??? task "`minimap2`: Read Alignment Details"' : '??? toggle "Map Reads to Host"'}) }}

{{ include_md("common_text/parse_mapping_task.md", condition="bam_to_unaligned_fastq", indent=4, replacements={'??? task "`parse_mapping`"' : '??? toggle "Extract Unaligned Reads"'}) }}

{{ include_md("common_text/assembly_metrics_task.md", indent=4, replacements={'??? task "`assembly_metrics`"' : '??? toggle "Host Read Mapping Statistics"'}) }}

    !!! techdetails "Host Decontaminate Technical Details"
        |  | Links |
        | --- | --- |
        | Subworkflow | [wf_host_decontaminate.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_host_decontaminate.wdl) |