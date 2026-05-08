---
title: Workflow Fragment `host_decontaminate`
fragment: true
---
<!-- if: theiaviral -->
??? task "`host_decontaminate`: Host Read Decontamination"

    Host genetic data is frequently incidentally sequenced alongside pathogens, which can negatively affect the quality of downstream analysis. Host Decontaminate attempts to remove host reads by aligning to a reference host genome that is directly inputted or acquired on-the-fly. The reference host genome can be inputted into the `host` input field as an assembly file (with `is_genome` set to "true"), acquired via [NCBI Taxonomy-compatible](https://www.ncbi.nlm.nih.gov/taxonomy) taxon input, or assembly accession (with `is_accession` set to "true"). Host Decontaminate maps inputted reads to the host genome using `minimap2`, reports mapping statistics to this host genome, and outputs the unaligned dehosted reads. 

    The detailed steps and tasks are as follows:

{{ include_md("common_text/estimate_genome_length_task.md", indent=4) }}

{{ include_md("common_text/ncbi_datasets_task.md", condition="theiaviral", indent=4) }}
<!-- endif -->

<!-- if: read_qc_trim -->
??? task "`read_decontaminate`: Mapping-based Read Decontamination"

    Known contaminant genetic data can be removed by mapping directly to an inputted `read_decontaminate_fasta`. This input can be a host genome, common microbial contaminant genome, or intentionally spiked sequences. The mapping statistics and aligned reads to the contaminant FASTA are outputted in JSON-formatted mappings, while downstream quality control tasks will input the decontaminated reads. An optional "pass/fail" status can be outputted based on identification of expected/unexpected sequences if the `expected_contaminants` input is populated with a comma-delimitted string of expected sequence headers - `expected_contaminants` must exactly match sequence headers in the input.

    The detailed steps and tasks are as follows:
<!-- endif -->

{{ include_md("common_text/minimap2_task.md", condition="only_map_ont", indent=4) }}

{{ include_md("common_text/parse_mapping_task.md", condition="bam_to_unaligned_fastq", indent=4, replacements={'??? task "`parse_mapping`: BAM File Handling"' : '??? task "`parse_mapping`: Extract Unaligned Reads"'}) }}

{{ include_md("common_text/mapping_stats_task.md", indent=4, replacements={'??? task "`mapping_stats`"' : '??? toggle "Host/Contaminant Read Mapping Statistics"'}) }}

{{ include_md("common_text/contaminant_check_task.md", indent=4, replacements={'??? task "`contaminant_check`"' : '??? toggle "Contaminant Detection Status"'}) }}

    !!! techdetails "Read Decontaminate Technical Details"
        |  | Links |
        | --- | --- |
        | Subworkflow | [wf_read_decontaminate.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_read_decontaminate.wdl) |
