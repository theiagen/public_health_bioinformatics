??? task "`megahit`"

<!-- if: theiaviral -->
    The `megahit` task is a wrapper for the [MEGAHIT assembler](https://github.com/voutcn/megahit), which is used for *de novo* metagenomic assembly of the cleaned reads. MEGAHIT is a fast and memory-efficient *de novo* assembler that can handle large datasets. This task is optional, turned off by default, and will only be called if MetaviralSPAdes fails. It can be enabled by setting the `skip_metaviralspades` parameter to `true`. The `megahit` task is used as a fallback option if the `spades` task fails during execution (see task `spades` for more details).

<!-- endif -->

    ???+ warning "Non-deterministic output(s)"
        This task may yield non-deterministic outputs.

    !!! techdetails "MEGAHIT Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_megahit.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_megahit.wdl) |
        | Software Source Code | [MEGAHIT on GitHub](https://github.com/voutcn/megahit) |
        | Software Documentation | [MEGAHIT](https://github.com/voutcn/megahit/blob/master/README.md) |
        | Original Publication(s) | [MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph](https://doi.org/10.1093/bioinformatics/btv033) |