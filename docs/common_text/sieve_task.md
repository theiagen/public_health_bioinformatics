---
title: Task Fragment `sieve`
fragment: true
---
??? task "`sieve`: Gene Screening and Typing"
    Sieve searches a sample's reads or assembly for a curated set of reference genes.

    The Sieve `engine` input defaults to BLAST for assemblies and KMA for raw reads, though KMA can also be used with assemblies. A gene hit is only kept if it matches the reference closely enough (`min_identity`) and over enough of the reference gene's length (`min_coverage`). Additional alignment settings can be passed as a comma-delimited list of `KEY=VALUE` pairs.

    Sieve reports a tab-delimited summary of the call, and can optionally output serogroup, a comma-delimited list of identified genes, and notes.

<!-- if: nmen -->
    **_Neisseria meningitidis_ serogrouping in TheiaProk**

    When GAMBIT identifies a sample as _Neisseria meningitidis_, TheiaProk automatically runs Sieve with the `nmen_serogroup` plugin (BLAST against the assembly) to predict the capsule serogroup.

    The assembly is screened against a database of _N. meningitidis_ capsule locus alleles, built from AllTheBacteria annotations of publicly available genomes and named according to the standard NEIS capsule locus nomenclature. Because capsule alleles are diverse and a capsule gene can be split across two contigs, the default thresholds are permissive (50% identity, 70% coverage), and hits from the same gene broken across a contig boundary are stitched back together before the gene is judged present or absent.

    The genes found are passed through a decision cascade ported from the PubMLST/BMGAP (PMGA) serogrouping logic. The cascade identifies which serogroup's capsule "backbone" is present, then confirms that every gene essential to that serogroup is present and intact — a gene that is fragmented, interrupted by an insertion element, or carrying a premature stop codon will prevent a positive call. The near-identical W and Y capsule polymerase genes (`csw` and `csy`) are resolved by translating the gene and reading the single amino acid that distinguishes W from Y.

    ??? toggle "Interpreting the `sieve_nmeningitidis_serogroup` output"
        | Value | Meaning |
        | --- | --- |
        | `A`, `B`, `C`, `W`, `X`, `Y`, `Z`, `E`, etc. | A confident serogroup call: the capsule backbone for that serogroup was found and all of its essential genes are present and intact. |
        | `NG` | Non-groupable. No capsule genes were found, an essential gene is missing, or the sample carries the capsule-null locus (`cnl`). |
        | `Inconclusive` | A capsule backbone was identified, but at least one essential gene is damaged (fragmented, disrupted, or truncated by an internal stop codon), so a serogroup cannot be confirmed. |
        | `Contaminated` | Capsule genes unique to more than one serogroup were detected, suggesting a mixed or contaminated sample. |

        The `sieve_nmeningitidis_notes` output explains the reasoning behind every call (for example, `"W backbone: All essential capsule genes intact and present"` or `"C backbone: ctrA fragmented (61.4% cov)"`), and `sieve_nmeningitidis_genes_present` lists the capsule genes that were recovered. When a call is `NG`, `Inconclusive`, or `Contaminated`, the notes should always be consulted before interpreting the result.

    The default thresholds can be adjusted with the optional `sieve_nmeningitidis_min_identity` and `sieve_nmeningitidis_min_coverage` inputs, and additional BLAST settings can be supplied with `sieve_nmeningitidis_parameters`.
<!-- endif -->

    !!! techdetails "sieve Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_sieve.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/task_sieve.wdl) |
        | Software Source Code | [sieve on GitHub](https://github.com/theiagen/sieve) |
        | Software Documentation | [sieve on GitHub](https://github.com/theiagen/sieve) |
        | Original Publication(s) | [Description and Nomenclature of _Neisseria meningitidis_ Capsule Locus](https://doi.org/10.3201/eid1904.111799) |
