---
title: Task Fragment `variant_annotate`
fragment: true
---
??? task "`variant_annotate`: Variant Effect Annotation"
    This task annotates variants in regions-of-interest using Ensembl Variant Effect Predictor (VEP) by reporting the predicted consequence of each variant (e.g. `missense_variant`, `frameshift_variant`) alongside its pseudo-[HGVS](https://hgvs-nomenclature.org/stable/) coding (`c.`) and protein (`p.`) notation.

    The task proceeds in three steps:

    1. **Variant extraction** - when `query_genes` and/or `query_genes_bed` are supplied, the variants overlapping the CDS coordinates of those genes are extracted into a sub-VCF, and each retained record is tagged with the overlapping gene name(s) in a `GENE` INFO field. Without either input, the full VCF is annotated.
    2. **Annotation** - VEP annotates the VCF against the reference FASTA and GFF by translating genes' CDS coordinates.
    3. **Reporting** - the VEP output is condensed into a comma-delimited list of human-readable annotations, where each transcript/protein identifier is replaced by the queried gene name, followed by the gene product name taken from the reference GFF.

    Each entry of the `variant_annotations` output takes the form:

    ```
    <query>: "<product>" (<consequence> <HGVSc> <HGVSp>; <ref>:<ref depth> <alt>:<alt depth>)
    ```

    For example:

    ```
    ERG11: "lanosterol 14-alpha demethylase" (missense_variant c.428A>G p.Lys143Arg; T:0 C:562)
    ```

    `<query>` is the `query_genes`/`query_genes_bed` term that selected the gene. When no query was supplied (the whole VCF is annotated) or no query matches the annotated feature, the product name is normalized into the label instead (`lanosterol.14-alpha.demethylase`). The product is quoted because product names frequently contain commas, which would otherwise be indistinguishable from the delimiter separating entries.

    Variants that resolve to neither a coding nor a protein change, and variants whose feature cannot be traced back to a gene product in the reference GFF, are omitted from `variant_annotations`. The full VEP run is preserved in `variant_annotation_summary` (VEP's HTML summary), while the extracted variants themselves are retained in `variant_annotation_gene_vcf`.

    ??? dna "Gene selection and coordinate sources"
        Variant annotation uses the same gene selection inputs as the `gene_coverage` task:

        - `query_genes` extracts gene coordinates from the reference GFF by matching the `product` qualifier of CDS entries. Matching is substring-based and case-insensitive unless `query_exact_match` is set to `true`
        - `query_genes_bed` supplies coordinates directly; gene names are taken from the fourth column of the BED file
        - If neither is supplied, the entire VCF is annotated

    ??? dna "Depth annotation with complementary nucleotides"
        Depths may be annotated with the complementary nucleotides that are associated with the single-nucleotide variant when the query region are annotated as the negative-strand.

    !!! techdetails "Variant Annotation Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_variant_annotate.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_variant_annotate.wdl) |
        | Software Source Code | [ensembl-vep on GitHub](https://github.com/Ensembl/ensembl-vep) |
        | Software Documentation | [VEP Documentation](https://useast.ensembl.org/info/docs/tools/vep/index.html) |
        | Original Publication(s) | [The Ensembl Variant Effect Predictor](https://doi.org/10.1186/s13059-016-0974-4) |
