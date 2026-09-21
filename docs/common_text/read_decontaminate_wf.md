---
title: Workflow Fragment `host_decontaminate`
fragment: true
---
<!-- if: theiaviral -->
??? task "`host_decontaminate`: Host Read Decontamination"

    Host genetic data is frequently incidentally sequenced alongside pathogens, which can negatively affect the quality of downstream analysis. Host Decontaminate attempts to remove host reads by aligning to a reference host genome that is directly input or acquired on-the-fly. The reference host genome can be provided into the `host` input field as an assembly file (with `is_genome` set to "true"), acquired via [NCBI Taxonomy-compatible](https://www.ncbi.nlm.nih.gov/taxonomy) taxon input, or assembly accession (with `is_accession` set to "true"). Host Decontaminate maps input reads to the host genome using `minimap2`, reports mapping statistics to this host genome, and outputs the unaligned dehosted reads.

    The detailed steps and tasks are as follows:

{{ include_md("common_text/estimate_genome_length_task.md", indent=4) }}

{{ include_md("common_text/ncbi_datasets_task.md", condition="theiaviral", indent=4) }}
<!-- endif -->

<!-- if: read_qc_trim -->
??? task "`mapped_read_removal`: Mapping-based Read Removal (optional)"
    Activate this task by providing a `mapped_read_removal_fasta`.

    Known contaminant genetic data can be removed by mapping directly to an provided `mapped_read_removal_fasta`. This input can be a host genome or a common microbial contaminant genome. The mapping statistics and aligned reads to the provided FASTA are created in JSON-formatted mappings, while downstream quality control tasks will input the reads that did not map. To additionally report a "pass/fail" status from expected/unexpected sequences, use `spike_in_screen`.

??? task "`spike_in_screen`: Mapping-based Spike-in Screening and Removal (optional)"
    Activate this task by providing **both** a `spike_in_fasta` and `expected_spike_ins`; the task will not run if either input is missing.

    Intentionally spiked sequences can be screened for by mapping directly to an provided `spike_in_fasta`. Reads that map to the `spike_in_fasta` are removed, and the remaining reads are created as `spike_in_removed_read1`/`spike_in_removed_read2` and passed to downstream quality control tasks. The mapping statistics and aligned reads to the spike-in FASTA are created in JSON-formatted mappings, alongside a "pass/fail" status based on identification of expected/unexpected sequences. `expected_spike_ins` is a comma-delimitted string of expected sequence headers that must exactly match sequence headers in the `spike_in_fasta`.

    If `mapped_read_removal` is also activated, `spike_in_screen` runs on the reads that `mapped_read_removal` left behind, so both sets of reads are removed before downstream analysis.

    The detailed steps and tasks are as follows:
<!-- endif -->

{{ include_md("common_text/minimap2_task.md", condition="only_map_ont", indent=4) }}

{{ include_md("common_text/parse_mapping_task.md", condition="bam_to_unaligned_fastq", indent=4, replacements={'??? task "`parse_mapping`: BAM File Handling"' : '??? task "`parse_mapping`: Extract Unaligned Reads"'}) }}

{{ include_md("common_text/mapping_stats_task.md", indent=4, replacements={'??? task "`mapping_stats`"' : '??? task "`mapping_stats`: Host/Contaminant Read Mapping Statistics"'}) }}

{{ include_md("common_text/contaminant_check_task.md", indent=4) }}

    !!! techdetails "Read Decontaminate Technical Details"
        |  | Links |
        | --- | --- |
        | Subworkflow | [wf_read_decontaminate.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_read_decontaminate.wdl) |
