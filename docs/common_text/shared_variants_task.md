---
title: Task Fragment `shared_variants`
fragment: true
---
<!-- if: snippy -->
??? task "`shared_variants`: Identifying Shared Variants (optional)"
    ##### Shared Variants (optional)

    ==This task activates when `call_shared_variants` is true.==
<!-- endif -->
<!-- if: find_shared_variants -->
??? task "`shared_variants`: Identifying Shared Variants"
    ##### Shared Variants
<!-- endif -->

    The Shared Variants task takes in the `concatenated_variants` output from the Concatenate Variants task and reshapes the data so that variants are rows and samples are columns. For each variant, samples where the variant was detected are populated with a "1" and samples where **either the variant was not detected or there was insufficient coverage to call variants** are populated with a "0". The resulting table is available as the `shared_variants_table` output.

    The `shared_variants_table` file will be in the following format:

    | CHROM | POS | TYPE | REF | ALT | FTYPE | STRAND | NT_POS | AA_POS | EFFECT | LOCUS_TAG | GENE | PRODUCT | sample1 | sample2 | sample3 |
    | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
    | PEKT02000007 | 2693938 | snp | T | C | CDS | - | 1008/3000 | 336/999 | synonymous_variant c.1008A>G p.Lys336Lys | B9J08_003879 | NA | chitin synthase 1 | 1 | 1 | 0 |
    | PEKT02000007 | 2529234 | snp | G | C | CDS | + | 282/336 | 94/111 | missense_variant c.282G>C p.Lys94Asn | B9J08_003804 | NA | cytochrome c | 1 | 1 | 1 |
    | PEKT02000002 | 1043926 | snp | A | G | CDS | - | 542/1464 | 181/487 | missense_variant c.542T>C p.Ile181Thr | B9J08_000976 | NA | dihydrolipoyl dehydrogenase | 1 | 1 | 0 |
    
    !!! techdetails "Shared Variants Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_shared_variants.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/utilities/task_shared_variants.wdl) |
