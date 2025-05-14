??? task "`bcftools_consensus`"

<!-- if: theiaviral -->
    The `bcftools_consensus` task generates a consensus genome by applying variants from the `clair3` task to a masked reference genome. It uses bcftools to filter variants based on the `min_depth` and `min_allele_freq` input parameter, left align and normalize indels, index the VCF file, and generate a consensus genome in FASTA format. Reference bases are substituted with filtered variants where applicable, preserved in regions without variant calls, and replaced with "N"s in areas masked by the `mask_low_coverage` task.

    ??? dna "`min_depth`"
        This parameter accepts an integer value to set the minimum read depth for variant calling and subsequent consensus sequence generation. The default value is `10`.

    ??? dna "`min_allele_freq`"
        This parameter accepts a float value to set the minimum allele frequency for variant calling and subsequent consensus sequence generation. The default value is `0.6`.
<!-- endif -->

    !!! techdetails "`bcftools_consensus` Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_bcftools_consensus.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_bcftools_consensus.wdl) |
        | Software Source Code | [bcftools on GitHub](https://github.com/samtools/bcftools) |
        | Software Documentation | [bcftools Manual Page](https://samtools.github.io/bcftools/bcftools.html) |
        | Original Publication(s) | [Twelve Years of SAMtools and BCFtools](https://doi.org/10.1093/gigascience/giab008) |