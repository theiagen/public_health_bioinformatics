---
title: Task Fragment `gatk_variants`
fragment: true
---
??? task "`gatk_variants`: Variant Calling (Illumina)"
    This task calls variants from an Illumina alignment (BAM) against the resolved reference genome using [GATK](https://gatk.broadinstitute.org/). 

    Variant calling consists of two GATK steps:

    - `HaplotypeCaller` performs local de-novo assembly of haplotypes in active (variable) regions to identify single-nucleotide polymorphisms (SNPs) and indels, emitting an intermediate per-sample GVCF.
    - `GenotypeGVCFs` converts genoetype likelihoods produced by the previous step into the final genotyped GVCF used for filtering.

    ??? dna "`gatk_ploidy` input parameter"
        Sample ploidy (N) passed to `HaplotypeCaller`. The default is `1` (haploid), which is appropriate for most default TheiaEuk use-cases. However, results can suffer if it is incorrectly called. Rerunning may be appropriate if a diploid or polyploid population is assembled via TheiaEuk.

    ??? dna "`gatk_intervals_file` input parameter"
        Optional intervals file (e.g. a BED or Picard-style interval list) restricting both `HaplotypeCaller` and `GenotypeGVCFs` to specific genomic regions. When omitted, the whole genome is analyzed.

    !!! techdetails "GATK Variants Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_gatk_variants.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_gatk_variants.wdl) |
        | Software Source Code | [GATK on GitHub](https://github.com/broadinstitute/gatk) |
        | Software Documentation | [GATK Documentation](https://gatk.broadinstitute.org/hc/en-us) |
        | Original Publication(s) | [A framework for variation discovery and genotyping using next-generation DNA sequencing data](https://doi.org/10.1038/ng.806) |
