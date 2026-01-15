<!-- if: theiaviral -->
??? task "`spades`"
<!-- endif -->
<!-- if: theiaprok -->
??? task "`SPAdes`: _De novo_ Assembly (alternative)"
    To activate this task, set `assembler` to `spades`.
<!-- endif -->

    `SPAdes` (St. Petersburg genome assembler) is a _de novo_ assembly tool that uses de Bruijn graphs to assemble genomes from Illumina short reads.

<!-- if: theiaprok -->
    In TheiaProk, SPAdes is run in `--isolate` mode, which is the recommended flag for high-coverage isolate and multi-cell Illumina data, which is typical of most bacterial sequencing projects. This method is optimized for improving assembly quality and decreasing runtime.
<!-- endif -->

<!-- if: theiaviral -->
    It is run with the `--metaviral` option, which is recommended for viral genomes. MetaviralSPAdes pipeline consists of three independent steps, `ViralAssembly` for finding putative viral subgraphs in a metagenomic assembly graph and generating contigs in these graphs, `ViralVerify` for checking whether the resulting contigs have viral origin and `ViralComplete` for checking whether these contigs represent complete viral genomes. For more details, please see the original publication.

    MetaviralSPAdes was selected as the default assembler because it produces the most complete viral genomes within TheiaViral, determined by CheckV quality assessment (see task `checkv` for technical details).

    ??? dna "`call_metaviralspades` input parameter"
        This parameter controls whether or not the `spades` task is called by the workflow. By default, `call_metaviralspades` is set to `true` because MetaviralSPAdes is used as the primary assembler. MetaviralSPAdes is generally recommended for most users, but it might not perform optimally on all datasets. If users encounter issues with MetaviralSPAdes, they can set the `call_metaviralspades` variable to `false` to bypass the `spades` task and instead *de novo* assemble using [MEGAHIT](https://github.com/voutcn/megahit) (see task `megahit` for details). Additionally, if the `spades` task fails during execution, the workflow will automatically fall back to using MEGAHIT for *de novo* assembly.
<!-- endif -->

    ???+ warning "Non-deterministic output(s)"
        This task may yield non-deterministic outputs.

    !!! techdetails "MetaviralSPAdes Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_spades.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_spades.wdl) |
        | Software Source Code | [SPAdes on GitHub](https://github.com/ablab/spades) |
        | Software Documentation | [SPAdes Manual](https://ablab.github.io/spades/index.html) |
        | Original Publication(s) | _TheiaProk_: [SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing](https://doi.org/10.1089/cmb.2012.0021)<br>_TheiaViral_: [MetaviralSPAdes: assembly of viruses from metagenomic data](https://doi.org/10.1093/bioinformatics/btaa490) |
