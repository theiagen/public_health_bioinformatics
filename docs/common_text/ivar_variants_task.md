??? task "`ivar_variants`"

    The `ivar_variants` task wraps the [iVar](https://andersen-lab.github.io/ivar/html/index.html) tool to call variants from the sorted BAM file produced by the `bwa` task. It uses the `ivar variants` command to identify and report variants based on the aligned reads. The `ivar_variants` task will filter all variant calls based on user-defined parameters, including `min_map_quality`, `min_depth`, and `min_allele_freq`. This task will return a VCF file containing the variant calls, along with the total number of variants, and the proportion of intermediate variant calls.

    ??? dna "`min_depth` input parameter"
        This parameter accepts an integer value to set the minimum read depth for variant calling and subsequent consensus sequence generation. The default value is `10`.

    ??? dna "`min_map_quality` input parameter"
        This parameter accepts an integer value to set the minimum mapping quality for variant calling and subsequent consensus sequence generation. The default value is `20`.

    ??? dna "`min_allele_freq` input parameter"
        This parameter accepts a float value to set the minimum allele frequency for variant calling and subsequent consensus sequence generation. The default value is `0.6`.

    !!! techdetails "iVar Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_ivar_variant_call.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_ivar_variant_call.wdl) |
        | Software Source Code | [Ivar on GitHub](https://andersen-lab.github.io/ivar/html/) |
        | Software Documentation | [Ivar Documentation](https://andersen-lab.github.io/ivar/html/manualpage.html) |
        | Original Publication(s) | [An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar](http://dx.doi.org/10.1186/s13059-018-1618-7) |