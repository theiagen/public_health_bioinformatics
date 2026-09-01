---
title: Task Fragment `gatk_filter`
fragment: true
---
??? task "`gatk_filter`: Variant Filtering (Illumina)"
    This task hard-filters the genotyped GVCF produced by the `gatk_variants` task using [GATK](https://gatk.broadinstitute.org/) `VariantFiltration`, then extracts the filtered records with `SelectVariants`.

    `VariantFiltration` does not remove variants; it annotates the `FILTER` column of each record with the name of any filter it fails (records passing all filters are marked `PASS`). Optional threshold inputs, depicted below, are automatically converted to filter expressions, though users can also input their own JEXL-formatted expression via the `gatk_filter_expression` input:

    | Input | Filter name | Expression | Flags a variant when… |
    | --- | --- | --- | --- |
    | `gatk_filter_min_variant_quality` | `variant_quality_filter` | `QUAL < value` | variant quality is below the threshold |
    | `gatk_filter_min_depth` | `depth_filter` | `DP < value` | read depth is below the threshold |
    | `gatk_filter_min_map_quality` | `mapping_quality_filter` | `MQ < value` | RMS mapping quality is below the threshold |
    | `gatk_filter_min_quality_by_depth` | `quality_by_depth_filter` | `QD < value` | quality-by-depth is below the threshold |

    Any threshold left unset is simply not applied. If none of the above are provided, `VariantFiltration` runs without threshold filters.

    ??? dna "`gatk_filter_expression` input parameter"
        A free-form [JEXL expression](https://gatk.broadinstitute.org/hc/en-us/articles/360035891011-JEXL-filtering-expressions) for custom filtering (annotated with the filter name `user_filter`. Use this for compound or non-standard criteria not covered by the named threshold inputs).

    !!! info "Outputs"
        - `gatk_filtered_vcf` — the GVCF with `FILTER`-column annotations applied by `VariantFiltration`.
        - `gatk_selected_vcf` — the records extracted by `SelectVariants`.

    !!! techdetails "GATK Filter Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_gatk_filter.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_gatk_filter.wdl) |
        | Software Source Code | [GATK on GitHub](https://github.com/broadinstitute/gatk) |
        | Software Documentation | [GATK VariantFiltration](https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration) |
        | Original Publication(s) | [A framework for variation discovery and genotyping using next-generation DNA sequencing data](https://doi.org/10.1038/ng.806) |
