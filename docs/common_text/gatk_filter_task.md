---
title: Task Fragment `gatk_filter`
fragment: true
---
??? task "`gatk_filter`: Variant Filtering (Illumina)"
    This task hard-filters the genotyped GVCF using [GATK](https://gatk.broadinstitute.org/) `VariantFiltration`, then extracts the filtered records with `SelectVariants`.

    !!! info "Outputs"
        - `gatk_filtered_vcf` — the GVCF with `FILTER`-column annotations applied by `VariantFiltration`. All records are reported.
        - `gatk_selected_vcf` — the records annotated with `PASS` extracted by `SelectVariants`.

    `VariantFiltration` annotates the `FILTER` column of each record with the name of any filter it fails (records passing all filters are marked `PASS`). Optional threshold inputs, depicted below, are automatically converted to filter expressions, though users can also input their own [JEXL-formatted expression](https://gatk.broadinstitute.org/hc/en-us/articles/360035891011-JEXL-filtering-expressions) via the `gatk_filter_expression` input:

    | Input | Filter name | Expression | Flags a variant when… |
    | --- | --- | --- | --- |
    | `gatk_filter_min_variant_quality` | `variant_quality_filter` | `QUAL < value` | variant quality is below the threshold |
    | `gatk_filter_min_depth` | `depth_filter` | `DP < value` | read depth is below the threshold |
    | `gatk_filter_min_map_quality` | `mapping_quality_filter` | `MQ < value` | RMS mapping quality is below the threshold |
    | `gatk_filter_min_quality_by_depth` | `quality_by_depth_filter` | `QD < value` | quality-by-depth is below the threshold |
    | `gatk_filter_max_fisher_strand_bias` | `fisher_strand_bias_filter` | `FS > value` | Fisher's strand bias exceeds the threshold |
    | `gatk_filter_max_strand_odds_ratio` | `strand_odds_ratio_filter` | `SOR > value` | the strand odds ratio exceeds the threshold |

    ??? dna "`gatk_filter_expression` input parameter"
        A free-form [JEXL expression](https://gatk.broadinstitute.org/hc/en-us/articles/360035891011-JEXL-filtering-expressions) for custom filtering (annotated with the filter name `user_filter`. Use this for compound or non-standard criteria not covered by the named threshold inputs).

    !!! techdetails "GATK Filter Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_gatk_filter.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_gatk_filter.wdl) |
        | Software Source Code | [GATK on GitHub](https://github.com/broadinstitute/gatk) |
        | Software Documentation | [GATK VariantFiltration](https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration) |
        | Original Publication(s) | [A framework for variation discovery and genotyping using next-generation DNA sequencing data](https://doi.org/10.1038/ng.806) |
