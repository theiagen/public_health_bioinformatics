---
title: Task Fragment `concatenate_variants`
fragment: true
---
<!-- if: snippy -->
??? task "`concatenate_variants`: Variant Concatenation (optional)"
    ##### Concatenate Variants (optional)

    ==This task activates when `call_shared_variants` is true.==
<!-- endif -->
<!-- if: find_shared_variants -->
??? task "`concatenate_variants`: Variant Concatenation"
    ##### Concatenate Variants
<!-- endif -->

    This task concatenates variant data from multiple samples into a single file `concatenated_variants`. It is very similar to the cat_files task, but also adds a column to the output file that indicates the sample associated with each row of data.

    The `concatenated_variants` file will be in the following format:

    | samplename | CHROM | POS | TYPE | REF | ALT | EVIDENCE | FTYPE | STRAND | NT_POS | AA_POS | EFFECT | LOCUS_TAG | GENE | PRODUCT |
    | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
    | sample1 | PEKT02000007 | 5224 | snp | C | G | G:21 C:0 |  |  |  |  |  |  |  |  |
    | sample2 | PEKT02000007 | 34112 | snp | C | G | G:32 C:0 | CDS | + | 153/1620 | 51/539 | missense_variant c.153C>G p.His51Gln | B9J08_002604 | hypothetical protein |  |
    | sample3 | PEKT02000007 | 34487 | snp | T | A | A:41 T:0 | CDS | + | 528/1620 | 176/539 | missense_variant c.528T>A p.Asn176Lys | B9J08_002604 | hypothetical protein |  |

    !!! techdetails "Concatenate Variants Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_cat_files.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/file_handling/task_cat_files.wdl) |
