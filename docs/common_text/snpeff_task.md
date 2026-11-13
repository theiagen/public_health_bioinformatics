---
title: Task Fragment `snpeff_gene_variants`
fragment: true
---
??? task "`snpeff_gene_variants`: Variant Effect Annotation"
    This task annotates variants in regions-of-interest using SnpEff. The predicted consequence of each variant is reported (e.g. `missense_variant`, `frameshift_variant`) alongside its [HGVS](https://hgvs-nomenclature.org/stable/)-like abbreviation.

    The task proceeds in four steps:

    1. **Database preparation** - reference GFF annotations and the reference genome FASTA are used to generate a custom SnpEff database.     
    2. **Variant extraction** - variants that overlap the CDS coordinates of `query_genes`/`query_genes_bed` are extracted and annotated into a sub-VCF.
    3. **Annotation** - SnpEff annotates the consequences of gene variants from the sub-VCF.
    4. **Reporting** - the SnpEff VCF annotations entries are transformed into a gene-labelled report

    Each entry of the `snpeff_gene_variants_hits` output takes the form:

    ```
    <query>: "<product>" (<consequence> <HGVSc> <HGVSp>; <ref>:<ref depth> <alt>:<alt depth>)
    ```

    For example:

    ```
    ERG11: "lanosterol 14-alpha demethylase" (missense_variant c.428A>G p.Lys143Arg; T:0 C:562)
    ```

    `<query>` is the `query_genes`/`query_genes_bed` term.

    The `snpeff_gene_variants_hits_abbreviated` output condenses each entry into `<query>:<amino acid change>` using one-letter residue codes (e.g. `FKS1:S639F`), falling back to `<query>:<nucleotide change>` (e.g. `FKKS1:428A>G`) for a variant with no protein change. When no variants are reported, both outputs read `No variants detected: <query_genes>`.

    Entries are omitted if: any of their consequence terms are listed in `suppressed_variants` (see below), they are not annotated against a transcript (e.g. `intergenic_region`, incomplete transcript, or lacking a stop codon), resolve to neither a coding nor a protein change, and if the variant cannot be traced back to a gene product in the reference GFF. The full set of annotations is retained in `snpeff_gene_variants_vcf`, alongside SnpEff's HTML summary (`snpeff_gene_variants_summary_html`). Every reported entry is also tabulated in `snpeff_gene_variants_report_tsv` as follows:

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

    ??? dna "`suppressed_variants` input parameter"
        `suppressed_variants` is a comma-delimited list of [SnpEff consequence terms](https://pcingola.github.io/SnpEff/snpeff/inputoutput/) (e.g. `upstream_gene_variant,downstream_gene_variant,synonymous_variant`) whose annotations are excluded from `snpeff_gene_variants_hits`, `snpeff_gene_variants_hits_abbreviated`, and `snpeff_gene_variants_report_tsv`. An annotation is excluded if any of its consequence terms are listed. By default, `upstream_gene_variant`, `downstream_gene_variant`, and `synonymous_variant` are suppressed. Suppressed annotations are still retained in `snpeff_gene_variants_vcf`.

        Supplying `suppressed_variants` replaces the defaults rather than adding to them, so re-list any defaults (e.g. `synonymous_variant`) that should still be excluded.

    ??? dna "Codon translation tables"
        By default, each contig is assigned the codon table from the GFF `transl_table` attribute. Contigs without a translation use SnpEff's default (Standard) table, and contigs with multiple translation tables within a single contig raise an error. A genome-wide translation table can be input via the `translation_table` input, referencing [NCBI's translation table integer codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).

    ??? dna "Depth annotation with complementary nucleotides"
        Depths are annotated with the reference and alternate nucleotides as they appear in the VCF (forward strand). When the query region is annotated on the negative strand, these are the complements of the nucleotides in the HGVS coding notation (e.g. `c.428A>G` with depths `T:0 C:562`).

    !!! techdetails "SnpEff Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_snpeff.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snpeff.wdl) |
        | Software Source Code | [SnpEff on GitHub](https://github.com/pcingola/SnpEff); [theiagene on GitHub](https://github.com/theiagen/theiagene) |
        | Software Documentation | [SnpEff Documentation](https://pcingola.github.io/SnpEff/) |
        | Original Publication(s) | [A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff](https://doi.org/10.4161/fly.19695) |
