??? task "`digger_denovo`: _De novo_ Assembly"
    _De novo_  assembly is the process or product of attempting to reconstruct a genome from scratch (without prior knowledge of the genome) using sequence reads. Assembly of fungal genomes from short-reads will produce multiple contigs per chromosome rather than a single contiguous sequence for each chromosome.

    In TheiaProk and TheiaEuk Illumina workflows, _de novo_ assembly is performed for samples that have sufficient read quantity and quality using [digger_denovo](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_digger_denovo.wdl), a subworkflow based off of [Shovill](https://github.com/tseemann/shovill) pipeline. The name "digger" is a nod to Shovill and SPAdes.
    
    ??? toggle "_De novo_ Assembly"
        !!! dna "`assembler` with `skesa` (default), `spades`, or `megahit`"
            To activate a particular assembler, set the `assembler` input parameter to either `skesa` (default), `spades`, or `megahit`.
        
            These tasks are mutually exclusive.

{{ include_md("common_text/skesa_task.md", indent=8) }}
{{ include_md("common_text/spades_task.md", indent=8, condition="theiaprok") }}
{{ include_md("common_text/megahit_task.md", indent=8, condition="theiaprok") }}

    ??? toggle "Assembly Polishing (optional)"
        To activate assembly polishing, set `call_pilon` to `true`.

{{ include_md("common_text/bwa_task.md", indent=8, condition="digger") }}
{{ include_md("common_text/pilon_task.md", indent=8, condition="digger") }}

    ??? toggle "Contig Filtering (optional)"

{{ include_md("common_text/filter_contigs_task.md", indent=8, condition="digger") }}

    !!! techdetails "Digger-Denovo Technical Details"
        |  | Links |
        | --- | --- |
        | Subworkflow | [wf_digger_denovo.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_digger_denovo.wdl) |
