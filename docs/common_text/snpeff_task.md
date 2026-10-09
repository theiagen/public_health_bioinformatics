---
title: Task Fragment `snpeff_gene_variants`
fragment: true
---
??? task "`snpeff_gene_variants`: Variant Effect Annotation"
    This task annotates variants in regions-of-interest using SnpEff by reporting the predicted consequence of each variant (e.g. `missense_variant`, `frameshift_variant`) alongside its pseudo-[HGVS](https://hgvs-nomenclature.org/stable/) coding (`c.`) and protein (`p.`) notation. A custom SnpEff database is built for each sample from the reference FASTA and GFF, so no pre-built SnpEff database is required.

    The task proceeds in four steps:

    1. **Database preparation** - the annotation section of the reference GFF (dropping any embedded `##FASTA` section) and the reference genome FASTA are staged as a custom SnpEff database. Each contig is assigned the codon table its GFF features declare via the `transl_table` attribute; contigs declaring none fall back to SnpEff's default (Standard) table.
    2. **Variant extraction** - when `query_genes` and/or `query_genes_bed` are supplied, the variants overlapping the CDS coordinates of those genes are extracted into a sub-VCF, and each retained record is tagged with the overlapping gene name(s) in a `GENE` INFO field. Without either input, the full VCF is annotated.
    3. **Annotation** - SnpEff builds the database and annotates the VCF against it. Upstream/downstream annotations are disabled, so only variants overlapping annotated features are annotated. When the VCF holds no variants, the build and annotation are skipped.
    4. **Reporting** - the SnpEff `ANN` entries are rendered into a gene-labelled report, where each entry is resolved back to its CDS product name through the reference GFF feature hierarchy.

    Each entry of the `snpeff_gene_variants_report` output takes the form:

    ```
    <query>: "<product>" (<consequence> <HGVSc> <HGVSp>; <ref>:<ref depth> <alt>:<alt depth>)
    ```

    For example:

    ```
    ERG11: "lanosterol 14-alpha demethylase" (missense_variant c.428A>G p.Lys143Arg; T:0 C:562)
    ```

    `<query>` is the `query_genes`/`query_genes_bed` term that selected the gene. When no query was supplied (the whole VCF is annotated) or no query matches the annotated feature, the product name is normalized into the label instead (`lanosterol.14-alpha.demethylase`). The product is quoted because product names frequently contain commas, which would otherwise be indistinguishable from the delimiter separating entries. The trailing per-allele read depths are taken from the VCF's `AD` (else `RO`/`AO`) field and are omitted when the record carries neither.

    The `snpeff_gene_variants_report_abbreviated` output condenses each entry into `<query>:<amino acid change>` using one-letter residue codes (e.g. `ERG11:K143R`), falling back to `<query>:<nucleotide change>` (e.g. `ERG11:428A>G`) for a variant with no protein change. When no variants are reported, both outputs read `No variants detected: <query_genes>`.

    Entries that are not annotated against a transcript (e.g. `intergenic_region`), entries that resolve to neither a coding nor a protein change, and entries whose variant cannot be traced back to a gene product in the reference GFF are omitted from the report. A transcript SnpEff flags as incomplete or lacking a start codon is treated as noncoding, so its protein change is not reported. The full set of annotations is retained in `snpeff_gene_variants_vcf`, alongside SnpEff's HTML summary (`snpeff_gene_variants_summary_html`) and per-gene effect counts (`snpeff_gene_variants_genes_txt`). Every reported entry is also tabulated in `snpeff_gene_variants_report_tsv`:

    | Column | Content | Example |
    | --- | --- | --- |
    | `GENE` | gene label | `ERG11` |
    | `HGVSc` | SnpEff's HGVS.c, prefixed with the transcript ID | `rna-x:c.428A>G` |
    | `HGVSp` | SnpEff's HGVS.p, prefixed with the CDS `protein_id` (else the transcript ID) | `prot-x:p.Lys143Arg` |
    | `NT` | abbreviated nucleotide change | `428A>G` |
    | `AA` | abbreviated amino acid change using one-letter residue codes | `K143R` |
    | `REPORT` | formatted report entry (see above) | `ERG11: "lanosterol 14-alpha demethylase" (...)` |

    Fields with no HGVS string are reported as `NA`.

    ??? dna "Gene selection and coordinate sources"
        Variant annotation uses the same gene selection inputs as the `gene_coverage` task:

        - `query_genes` extracts gene coordinates from the reference GFF by matching the `product` or `locus_tag` qualifier of gene/CDS entries. Matching is substring-based and case-insensitive unless `query_exact_match` is set to `true`
        - `query_genes_bed` supplies gene names from the fourth column of the BED file when `query_genes` is not supplied; coordinates are always taken from the reference GFF
        - If neither is supplied, the entire VCF is annotated

    ??? dna "Depth annotation with complementary nucleotides"
        Depths are annotated with the reference and alternate nucleotides as they appear in the VCF (forward strand). When the query region is annotated on the negative strand, these are the complements of the nucleotides in the HGVS coding notation (e.g. `c.428A>G` with depths `T:0 C:562`).

    ??? dna "Deviations from HGVS"
        HGVS strings are SnpEff's, so reported changes follow the [HGVS recommendations](https://hgvs-nomenclature.org/stable/) except where noted below. One-letter amino acid codes, `*` for a stop codon, and the short frameshift form (`p.Lys128fs`) are permitted by HGVS and are not deviations.

        | Deviation | HGVS | Reported |
        | --- | --- | --- |
        | reference sequence identifier dropped (the gene label stands in for it) | `NM_000001.1:c.428A>G` | `c.428A>G` |
        | parentheses around predicted protein changes absent | `p.(Lys143Arg)` | `p.Lys143Arg` |
        | synonymous change repeats the reference residue instead of using `=` | `p.Asp164=` | `p.Asp164Asp`, `D164D` |
        | deleted/duplicated bases listed | `c.383del` | `c.383delA` |
        | protein-level duplication described as an insertion | `p.Gly294_Ser297dup` | `p.Ala292_Ser293insSerGlySerAla` |
        | protein-level indel not shifted 3′ through a repeat | `p.Gln444_Gly448del` | `p.Phe432_Gly436del` |

    !!! techdetails "SnpEff Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_snpeff.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snpeff.wdl) |
        | Software Source Code | [SnpEff on GitHub](https://github.com/pcingola/SnpEff); [theiagene on GitHub](https://github.com/theiagen/theiagene) |
        | Software Documentation | [SnpEff Documentation](https://pcingola.github.io/SnpEff/) |
        | Original Publication(s) | [A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff](https://doi.org/10.4161/fly.19695) |
