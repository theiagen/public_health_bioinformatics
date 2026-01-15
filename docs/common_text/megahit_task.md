<!-- if: theiaviral -->
??? task "`megahit`"
<!-- endif -->
<!-- if: theiaprok -->
??? task "`MEGAHIT`: _De novo_ Assembly (alternative)"
    To activate this task, set `assembler` to `megahit`.
<!-- endif -->

    The [MEGAHIT assembler](https://github.com/voutcn/megahit) is a fast and memory-efficient _de novo_ assembler that can handle large datasets. While optimized for metagenomics, MEGAHIT also performs well on single-genome assemblies, making it a versatile choice for various assembly tasks.

    MEGAHIT uses a multiple k-mer strategy that can be beneficial for assembling genomes with varying coverage levels, which is common in metagenomic samples. It constructs succinct de Bruijn graphs to efficiently represent the assembly process, allowing it to handle large and complex datasets with reduced memory usage.

<!-- if: theiaviral -->
    This task is optional, turned off by default, and will only be called if MetaviralSPAdes fails. It can be enabled by setting the `skip_metaviralspades` parameter to `true`. The `megahit` task is used as a fallback option if the `spades` task fails during execution (see task `spades` for more details).
<!-- endif -->

    ???+ warning "Non-deterministic output(s)"
        This task may yield non-deterministic outputs.

    !!! techdetails "MEGAHIT Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_megahit.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_megahit.wdl) |
        | Software Source Code | [MEGAHIT on GitHub](https://github.com/voutcn/megahit) |
        | Software Documentation | [MEGAHIT on GitHub](https://github.com/voutcn/megahit/blob/master/README.md) |
        | Original Publication(s) | [MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph](https://doi.org/10.1093/bioinformatics/btv033) |
