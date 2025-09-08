??? task "`ivar_variants`: Variant Calling"
    iVar uses the outputs of `samtools mpileup` to call single nucleotide variants (SNVs) and insertions/deletions (indels). Several key parameters can be set to determine the stringency of variant calling, including minimum quality, minimum allele frequency, and minimum depth.

    This task returns a VCF file containing all called variants, the number of detected variants, and the proportion of those variants with allele frequencies between 0.6 and 0.9 (also known as _intermediate_ variants).

<!-- if: theiacov -->
    For TheiaCoV, the following default parameters are used:
    
    - minimum quality: 20
    - minimum depth: 100
    - minimum allele frequency: 0.06
<!-- endif -->

<!-- if: theiaviral -->
    ??? dna "`min_depth` input parameter"
        This parameter accepts an integer value to set the minimum read depth for variant calling and subsequent consensus sequence generation. The default value is `10`.

    ??? dna "`min_map_quality` input parameter"
        This parameter accepts an integer value to set the minimum mapping quality for variant calling and subsequent consensus sequence generation. The default value is `20`.

    ??? dna "`min_allele_freq` input parameter"
        This parameter accepts a float value to set the minimum allele frequency for variant calling and subsequent consensus sequence generation. The default value is `0.6`.
<!-- endif -->
    !!! techdetails "iVar Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_ivar_variant_call.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_ivar_variant_call.wdl) |
        | Software Source Code | [Ivar on GitHub](https://andersen-lab.github.io/ivar/html/) |
        | Software Documentation | [Ivar Documentation](https://andersen-lab.github.io/ivar/html/manualpage.html) |
        | Original Publication(s) | [An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar](http://dx.doi.org/10.1186/s13059-018-1618-7) |