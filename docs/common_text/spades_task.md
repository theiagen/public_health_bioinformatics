??? task "`spades`"

<!-- if: theiaviral -->
    The `spades` task is a wrapper for the [SPAdes assembler](https://github.com/ablab/spades), which is used for *de novo* assembly of the cleaned reads. It is run with the `--metaviral` option, which is recommended for viral genomes. MetaviralSPAdes pipeline consists of three independent steps, `ViralAssembly` for finding putative viral subgraphs in a metagenomic assembly graph and generating contigs in these graphs, `ViralVerify` for checking whether the resulting contigs have viral origin and `ViralComplete` for checking whether these contigs represent complete viral genomes. For more details, please see the original publication.

    MetaviralSPAdes was selected as the default assembler because it produces the most complete viral genomes within TheiaViral, determined by CheckV quality assessment (see task `checkv` for technical details).

    ??? dna "`call_metaviralspades`"
        This parameter controls whether or not the `spades` task is called by the workflow. By default, `call_metaviralspades` is set to `true` because MetaviralSPAdes is used as the primary assembler. MetaviralSPAdes is generally recommended for most users, but it might not perform optimally on all datasets. If users encounter issues with MetaviralSPAdes, they can set the `call_metaviralspades` variable to `false` to bypass the `spades` task and instead *de novo* assemble using [MEGAHIT](https://github.com/voutcn/megahit) (see task `megahit` for details). Additionally, if the `spades` task fails during execution, the workflow will automatically fall back to using MEGAHIT for *de novo* assembly.

    ???+ warning "Important"
        In this workflow, *de novo* assembly is primarily used to facilitate the selection of a closely related reference genome, though high quality *de novo* assemblies can be used for downstream analysis. If the user provides an input `reference_fasta`, all subsequent assembly and reference selections tasks will be skipped, including:

        - `spades`
        - `megahit`
        - `checkv_denovo`
        - `quast_denovo`
        - `skani`
        - `ncbi_datasets`

    !!! techdetails "MetaviralSPAdes Technical Details"
    |  | Links |
    | --- | --- |
    | Task | [task_spades.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_spades.wdl) |
    | Software Source Code | [SPAdes on GitHub](https://github.com/ablab/spades) |
    | Software Documentation | [SPAdes Manual](https://ablab.github.io/spades/index.html) |
    | Original Publication(s) | [MetaviralSPAdes: assembly of viruses from metagenomic data](https://doi.org/10.1093/bioinformatics/btaa490) |
<!-- endif -->