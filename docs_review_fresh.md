# Documentation Review — Fresh Analysis

_Comprehensive analysis of all documentation in `docs/` for typos/grammar, accessibility, and ease-of-use. Generated fresh from the current state of the files, independent of any prior review._

> **Refreshed 2026-07-08.** Since this review was generated, the glossary/tooltip system and documentation linting tooling have been added and grammar/correctness fixes applied (committed in `4c81cfd1 "add a glossary and fix various typos"`). See **[Refresh status](#refresh-status-2026-07-08)** below for what has been resolved; the original analysis is preserved beneath it.

## Summary

- **Total findings:** 678 across 214 files
- **By category:** 27 typos · 166 grammar · 105 accessibility · 380 ease-of-use
- **By severity:** 8 high · 174 medium · 496 low

_(The counts above are the original snapshot; see Refresh status for what has since been resolved.)_

## Refresh status (2026-07-08)

Since generation, two systems were added and a round of fixes applied:

- **Glossary + tooltips** — new `docs/guides/glossary.md` (7 sections: Sequencing and analysis terms, File formats, Genomic characterization terms, Platforms and tools, Databases/repositories/organizations, Pathogens and organisms, Acronyms and abbreviations) plus `includes/abbreviations.md` auto-appended site-wide for organization tooltips.
- **Linting tooling** — `.vale.ini` + `styles/PHB/` (tool-name consistency, American spelling), `.markdownlint.yaml`, `.pre-commit-config.yaml` (whitespace/`inputted`/space-before-%), and `.github/workflows/docs-lint.yml`.
- **Content fixes** — grammar pass plus most of the high-priority correctness items (below).

### Systemic patterns — current status

| # | Pattern | Status | Mechanism |
| --- | --- | --- | --- |
| 1 | Undefined jargon & unexpanded acronyms | ✅ Addressed | Glossary page + auto-appended org tooltips; "expand on first use" house rule in `doc_contribution.md` |
| 2 | Dense paragraphs → lists | ⬜ Outstanding | Content-level; style rule added to contributing guide |
| 3 | Missing "why/when to use" framing | 🟡 Partial | Some fixed (e.g. `phylovalidate_task.md`); broadly outstanding |
| 4 | Tables with terse headers / no legend | ⬜ Outstanding | Content-level |
| 5 | Inconsistent tool-name capitalization | ✅ Tooling | Vale `PHB.ToolNames` (existing instances need one cleanup pass) |
| 6 | "inputted" / passive phrasing | ✅ Tooling | pre-commit `no-inputted` hook |
| 7 | Whitespace & punctuation hygiene | ✅ Tooling | pre-commit trailing-whitespace / EOF / space-before-% + markdownlint |
| 8 | Heading-hierarchy skips | ✅ Tooling | markdownlint `MD001` |
| 9 | Vague link text / color-only refs | 🟡 Partial | markdownlint `MD034` (bare URLs); "click here"/color-only still manual |
| 10 | Inconsistent terminology | 🟡 Partial | Some reconciled; broadly outstanding |
| 11 | Duplicated content w/ same defect | ⬜ Outstanding | Content-level |
| 12 | British → American spelling | ✅ Tooling | Vale `PHB.AmericanSpelling` |

Tooling now covers patterns **5, 6, 7, 8, 12** (and partially 9); the glossary covers **1**. Patterns **2, 3, 4, 10, 11** remain content-level work not amenable to automation.

### Highest-priority individual fixes — current status

| # | Item | Status |
| --- | --- | --- |
| 1 | `read_screen_task.md` `max_genome_size` logic inversion | ✅ Fixed ("larger than") |
| 2 | `augur_clades_task.md` `\alt` → `\talt` | ✅ Fixed |
| 3 | `nextclade_task.md` rabies species name | ✅ Fixed (_Lyssavirus rabies_) |
| 4 | `bwa_task.md` "Burrow-Wheeler" → "Burrows-Wheeler" | ✅ Fixed |
| 5 | `augur_refine_task.md` "TimeTree" → "TreeTime" | ✅ Fixed |
| 6 | `terra_2_ncbi.md` (line 96) misplaced `**` bold | ✅ Fixed (now **nothing** will appear) |
| 7 | `pangolin_update.md` garbled sentence | ✅ Fixed (rewritten) |
| 8 | `standalone/kraken2.md` table | ✅ Fixed (dates now ISO `YYYY-MM-DD`; EuPathDB48 rows disambiguated by version) |
| 9 | `pangolin_task.md` "Reoccuring" | ✅ Fixed |
| 10 | Plain-language framing: `extract_clade_mutations_task.md` & `phylovalidate_task.md` | ✅ Fixed (both now open with a plain-language lead) |
| 11 | `data_summary_task.md` TRUE/FALSE inconsistency | ✅ Fixed |
| 12 | `freyja.md` stray `=` and mismatched quote | ✅ Fixed |

**All twelve high-priority individual fixes are now resolved.** Remaining work is the content-level systemic patterns (#2 dense paragraphs, #4 table legends, #11 duplicated content, and the partial items #3/#9/#10-class framing across other files).

> The detailed per-file findings below are the **original snapshot** and have not been individually re-verified. Mechanical items (whitespace, tool names, British spelling, heading skips) are now caught by the linting tooling, and jargon/acronym items are mitigated by the glossary; the remaining prose-quality items (dense paragraphs, table legends, framing) are still valid as a content backlog.

## Systemic Patterns & Priorities

# Systemic Analysis: PHB Documentation Review

This review surfaced ~450 findings across ~200 files. The overwhelming majority are instances of a small number of recurring patterns. Fixing these at the standard/template level (rather than file-by-file) will resolve the bulk of the issues and prevent recurrence.

## Cross-Cutting Patterns

### 1. Undefined jargon & unexpanded acronyms on first use
**Widespread — by far the dominant pattern (~150+ findings, present in nearly every `common_text` fragment and most workflow pages).**
The same core terms recur unexplained: `k-mer`, `contig`, `de novo`, `coverage`/`depth`, `in silico`, `BAM`, `VCF`, `FASTA`/`FASTQ`, `BED`, `GFF3`, `SNP`/`SNV`, `indel`, `clade`, `serotype`/`serogroup`/`biotype`, `locus`, `allele`, and acronyms `ONT`, `AMR`, `MLST`, `WGS`, `QC`, `PE`/`SE`, `HRRT`, `ANI`, `RefSeq`, `GTDB`.
Most affected: the entire `docs/common_text/*` set; notably `amrfinderplus_task.md`, `kraken2_task.md`, `read_screen_task.md`, `theiaviral_inputs.md`, `metabuli_task.md`, `kaptive_task.md`, plus `guides/phylogenetics.md`.
**Recommended solution:** Build a single shared glossary page and (a) link core terms to it, and (b) adopt a house rule that each acronym is expanded on first use per page. Because these fragments share a fixed vocabulary, a one-time glossary + a "define once" convention resolves the majority of the review. A per-fragment glossary macro/snippet would let definitions live in one place.

### 2. Dense single paragraphs that should be lists / missing plain-language summary
**Widespread (~40 findings).** Long run-on paragraphs packing purpose + parameters + defaults + edge cases into one block.
Most affected: `data_summary_task.md`, `bcftools_consensus_task.md`, `amrfinderplus_task.md`, `qc_check_task.md`, `read_screen_task.md`, `midas_task.md`, `fastp_task.md`, `theiaviral_inputs.md`, `pangolin_task.md`, `custom_organisms.md`, `freyja.md`.
**Recommended solution:** Adopt a template pattern for task/parameter descriptions: one plain-language "what this does / when to use it" lead sentence, then a bulleted list for steps, parameters, and defaults. Add this to `doc_contribution.md` as a style rule.

### 3. Missing "why/when to use this" framing
**Common (~20 findings).** Fragments describe mechanics but never the user-facing purpose.
Most affected: `contaminant_check_task.md`, `dnaapler_task.md`, `download_terra_table_task.md`, `ksnp_task.md`, `snp_dists_task.md`, `find_shared_variants.md`, `rename_fastq.md`, `gambit_query.md`, `cauris_cladetyper.md`, `vadr_update.md`, `tbprofiler_tngs.md`.
**Recommended solution:** Add a required "Purpose / When to use" opening line to the workflow/task page template.

### 4. Tables with abbreviated columns and no legend/lead-in
**Common (~15 findings), an accessibility issue.** Example data tables use terse headers (CHROM, POS, FTYPE, NT_POS, AA_POS, LOCUS_TAG, StDev, CDS) or glyph-only cells (check/cross/plus) with the key placed after the table or absent.
Most affected: `concatenate_variants_task.md`, `shared_variants_task.md`, `gene_coverage_task.md`, `midas_task.md`, `create_terra_table.md`, `mercury_prep_n_batch.md`, `freyja.md`, `theiacov.md`, `arln_stats_task.md`, `sra_fetch.md`.
**Recommended solution:** Standard rule: every example/data table gets a one-sentence lead-in and a legend defining columns/symbols _before_ the table. Provide a reusable column-glossary snippet for the recurring variant-table headers.

### 5. Inconsistent tool-name capitalization / spelling
**Common (~15 findings).** The same tool is spelled multiple ways within or across files.
Instances: `Gubbins`/`gubbins`, `RAxML`/`RaxML`, `PBPtyper`/`PBPTyper`, `StxTyper`/`Stxtyper`, `Skani`/`skani`, `Docker`/`docker`, `miniwdl`/`miniWDL`, `GitHub`/`Github`, `IQ-TREE`/`IQTREE`/`IQTree`, `ModelFinder`/`Model Finder`, `Kraken2 database`/`Kraken DB`, `BaseSpace`/`Basespace`.
Most affected: `gubbins_task.md`, `snippy_streamline.md`/`snippy_streamline_fasta.md`, `theiaviral.md`, `commandline.md`, `stxtyper_task.md`, `pbptyper_task.md`.
**Recommended solution:** Maintain a canonical tool-name list in the contributing guide and add a CI/linter dictionary check (or vale style) for approved spellings.

### 6. Nonstandard "inputted" and awkward passive input phrasing
**Common (~12 findings).** "inputted" / "previously inputted" / "directly inputted" recur.
Most affected: `augur_export_task.md`, `contaminant_check_task.md`, `cophylogeny_task.md`, `estimate_genome_length_task.md`, `ete4_identify_task.md`, `extract_clade_mutations_task.md`, `read_decontaminate_wf.md`.
**Recommended solution:** Global find-and-replace "inputted" → "input"/"provided"; add to linter word-blocklist.

### 7. Whitespace & punctuation hygiene (trailing spaces, double spaces, percent spacing, missing terminal periods, missing comma after intro clause, bullet-punctuation inconsistency, missing table pipes)
**Very widespread (~50 findings), mostly low severity.** Includes `95 %`/`25 %`/`84.35 %` (space before `%`) vs `95%`, missing serial commas, and inconsistent bullet-list terminal punctuation.
Most affected: `gubbins_task.md`, `kraken2_task.md`, `kaptive_task.md`, `bbduk_task.md`, `theiaviral_inputs.md`, `centroid_task.md`, `standalone/kraken2.md`, `workflows_overview/*`.
**Recommended solution:** Add automated tooling — markdownlint + a trailing-whitespace/double-space pre-commit hook + a style rule "no space before `%`, use Oxford comma." This mechanically clears a large low-severity tail.

### 8. Heading hierarchy skips & heading style (accessibility)
**Several files.** H1→H3 jumps and full-sentence headings hinder screen-reader navigation.
Most affected: `workflows_overview/workflows_kingdom.md`, `workflows_overview/workflows_type.md`, `getting_started/commandline.md`, `guides/phylogenetics.md`.
**Recommended solution:** Enforce no-skipped-heading-levels via markdownlint (MD001) and template headings as noun phrases.

### 9. Vague link text, bare URLs, and color/emoji-only UI references (accessibility)
**Common (~12 findings).** Links like "this file", "this output format.", "please see this", bare URLs, and instructions that reference "the orange arrow", "✅ box", "two squares", or 📺 emoji.
Most affected: `stxtyper_task.md`, `data_summary_task.md`, `tbp_parser_task.md`, `mercury_prep_n_batch.md`, `guides/gambit.md`, `rasusa.md`, `ksnp3.md`/`ksnp4.md`, `ont_barcode_concatenation.md`.
**Recommended solution:** Style rule: descriptive link text (never "this file"/bare URL), and always pair color/icon references with a textual label (e.g., "the Execute (orange) button").

### 10. Inconsistent terminology for the same concept
**Several high-impact instances.** Same idea named two ways, which erodes trust and can mislead.
Instances: `min_genome_size`/`max_genome_size` vs `min_genome_length`/`max_genome_length` (`read_screen_task.md`); `SNP` vs `SNV` (`clair3_variants.md`, `snippy`); WDL = "Workflow Development" vs "Workflow Description" Language (`commandline.md`); `TimeTree` vs `TreeTime` (`augur_refine_task.md`); "Excel" vs "TSV" (`theiavalidate.md`); fungi vs "mycotics" (`index.md`).
**Recommended solution:** Pick one canonical term per concept (document in contributing guide) and reconcile. Prioritize `read_screen_task.md` variable names and the WDL expansion since these directly confuse users.

### 11. Duplicated / copy-pasted content carrying the same defect
**Several instances.** Shared text propagates errors.
Instances: the mislabeled "Resistance Gene Database" toggle appears identically in `abricate_vibrio_task.md` and `srst2_task.md` (both also have the "presence or absence...are used" agreement error); `ivar_consensus_task.md` parameter descriptions copied from `ivar_variants_task.md` (wrong task behavior); "Read Quantification" text reused in `fastqc_task.md` from `fastq_scan_task.md` mislabeling FastQC's purpose; `ksnp3.md`/`ksnp4.md` share identical "etc."/"else" grammar issues.
**Recommended solution:** Where fragments are near-duplicates, factor into a shared include; otherwise fix all copies together and note the linkage.

### 12. British→American spelling inconsistency
**Minor but recurring.** `characterisation` (`hicap_task.md`), `analyses` as a verb (`sonneityper_task.md`), `Orthologue` (`busco_task.md`).
**Recommended solution:** Add American-spelling dictionary to the linter (note genuine quotes/official names as exceptions).

## Highest-Priority Individual Fixes

These are high-severity one-offs (mostly correctness/factual) that don't fold into a pattern and should be fixed directly:

1. **`read_screen_task.md` (line 30) — logic inversion.** `max_genome_size` rationale says a sample fails if the estimated size is "smaller than" the maximum; it must be "larger than." A real correctness bug in documented behavior.
2. **`augur_clades_task.md` (line 6) — malformed format string.** Column separator `\alt` should be `\talt` (tab). Users copying the header will produce a broken `clades_tsv`.
3. **`nextclade_task.md` (line 11) — likely wrong species name.** `_L. rabies_` is non-standard; verify/correct the rabies virus taxonomic name.
4. **`bwa_task.md` (all occurrences) — factual name error.** "Burrow-Wheeler Aligner" should be "Burrows-Wheeler Aligner."
5. **`augur_refine_task.md` (line 6) — likely wrong tool name.** Body says "TimeTree" but the tech-details table and sibling Augur tasks use "TreeTime."
6. **`terra_2_ncbi.md` (line 96) — broken markup producing a garbled sentence.** Misplaced `**` bold markers ("...behave exactly **like a real submission, but since it's detached, nothing **will appear...").
7. **`pangolin_update.md` (line 9) — garbled core sentence** describing what the workflow does (docker-image lineage update); needs a full rewrite for comprehensibility.
8. **`standalone/kraken2.md` (lines 36–43) — misaligned table row + ambiguous dates.** The "viral w/human" row's columns are shifted (missing Source cell); dates use ambiguous D/M/Y (`12/1/2024`) inconsistent with ISO elsewhere; two rows share the name "EuPathDB48."
9. **`pangolin_task.md` (line 8) — misspelling in the scorpio expansion:** "Reoccuring" (verify correct scorpio wording).
10. **`extract_clade_mutations_task.md` & `phylovalidate_task.md` — no plain-language framing at all** for an explicitly non-technical audience (both rated high on ease-of-use); add a purpose sentence before the technical detail.
11. **`data_summary_task.md` — TRUE/empty vs FALSE inconsistency.** Prose says presence=TRUE/absence=empty, but the example CSV shows FALSE; reconcile so the example matches the description.
12. **`freyja.md` — copy/format artifacts:** stray trailing `=` after the `Freyja_Dashboard_PHB` link (line 34) and a mismatched quote in the `freyja_demixed` example (`0.1'`, line 158).

## Suggested Execution Order
1. Ship automated tooling first (markdownlint + spelling/word-blocklist + whitespace hook): clears patterns 5, 6, 7, 8, 12 and most of the low-severity volume mechanically.
2. Fix the 12 high-priority individual defects (correctness/factual).
3. Build the shared glossary + "expand acronyms on first use" convention (pattern 1) — the single largest comprehension win.
4. Update the page/task template with mandatory "purpose" line, bulleted parameter lists, and table legend/link-text rules (patterns 2, 3, 4, 9), then reconcile terminology and de-duplicated fragments (patterns 10, 11).

## All High-Severity Findings

| File | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| `docs/common_text/augur_clades_task.md` | typo | Line 6 | Malformed/escaped header string: 'header "clade\tgene\tsite\alt"' — the last separator is '\alt' (a single backslash) instead of '\t' (tab), inconsistent with the others. | Correct to '\talt' so the column separators are consistent (clade\tgene\tsite\talt). |
| `docs/common_text/data_summary_task.md` | ease-of-use | Line 12 | Very long, dense paragraph mixing input format, behavior, output, and an optional setting; difficult for a non-technical reader to follow. | Break into bullets: (1) what to enter in data_summary_* and sample_names, (2) input format example, (3) what the task produces (CSV), (4) the phandango_coloring option. |
| `docs/common_text/extract_clade_mutations_task.md` | ease-of-use | Line 6 | Highly technical with no plain-language framing: 'Augur-compatible clades.tsv', 'signature clade-defining sequences', 'nucleotide JSON outputted by Augur Ancestral', 'amino acid JSON outputted by Augur Translate' will be opaque to non-technical readers. | Add an opening sentence in plain language explaining the purpose (identifying the mutations that define each clade/group in a tree) before the technical details, and gloss Augur, JSON, clade. |
| `docs/common_text/nextclade_task.md` | typo | Line 11 | "_L. rabies_" appears to be an incorrect species abbreviation; rabies virus is _Lyssavirus rabies_ but "_L. rabies_" is non-standard/possibly wrong, and the linked dataset is for rabies. This may be a factual/typo error. | Verify and correct the species name (e.g., "rabies virus" or the correct italicized taxonomic name). |
| `docs/common_text/phylovalidate_task.md` | ease-of-use | Line 6 and Line 8 | The description is highly technical (polytomies, topologies, Lin-Rajan-Moret distance, matching cluster distance, Robinson-Foulds distance) with no plain-language explanation of what the tool does for the user or when they would use it. The audience is explicitly non-technical. | Add an opening plain-language sentence describing the purpose, e.g. "PhyloValidate checks whether two phylogenetic trees (evolutionary relationship diagrams) are essentially the same shape, which is useful for confirming that results are reproducible." Then keep the technical detail for advanced users. |
| `docs/common_text/read_screen_task.md` | ease-of-use | Line 30 (theiaviral table, max_genome_size row) | The rationale for `max_genome_size` reads 'A sample will fail the read screening if the estimated genome size is smaller than `max_genome_size`'. For a maximum threshold this is logically wrong — it should fail when the size is LARGER than the maximum. | Change 'smaller than' to 'larger than' for the max_genome_size row. |
| `docs/common_text/read_screen_task.md` | ease-of-use | Line 42 (theiacov table, max_genome_length row) / Line 54 (theiaprok) / Line 66 (theiaeuk) / Line 79 (theiaeukont) | Inconsistent variable naming between sections: the intro bullets and theiaviral table use `min_genome_size`/`max_genome_size`, while the theiacov/theiaprok/theiaeuk tables use `min_genome_length`/`max_genome_length`. Same concept, two different names, which is confusing. | Standardize on one term (e.g. `min_genome_length`/`max_genome_length`) across all sections, or note explicitly that the two names refer to the same input. |
| `docs/guides/phylogenetics.md` | ease-of-use | Line 45 ('heterologous sites'); also 'HGT' on Line 61 | Undefined jargon for a non-technical audience. 'heterologous sites' is highly technical and undefined; 'HGT' is an unexpanded acronym on first use. | Define or simplify 'heterologous sites' (e.g., 'sites where the sequence differs from the reference'). Expand HGT to 'horizontal gene transfer (HGT)' on first use. |

## All Typos

| File | Location | Issue | Suggestion |
| --- | --- | --- | --- |
| `docs/common_text/augur_clades_task.md` | Line 6 | Malformed/escaped header string: 'header "clade\tgene\tsite\alt"' — the last separator is '\alt' (a single backslash) instead of '\t' (tab), inconsistent with the others. | Correct to '\talt' so the column separators are consistent (clade\tgene\tsite\talt). |
| `docs/common_text/augur_refine_task.md` | Line 6 | Likely inconsistency: the surrounding Augur tasks describe time-calibration via 'TreeTime', but this file says 'TimeTree'. If the underlying software is TreeTime, this is a naming error. | Verify and correct the software name to 'TreeTime' if that is the tool used (the techdetails table on lines 12-14 references TreeTime, not TimeTree, confirming the body text is likely wrong). |
| `docs/common_text/bwa_task.md` | Lines 8, 11, 15, 19, 22 (all 'Burrow-Wheeler Aligner') | BWA is consistently expanded as 'Burrow-Wheeler Aligner'; the correct name is 'Burrows-Wheeler Aligner' (with an 's'). | Correct all occurrences to 'Burrows-Wheeler Aligner' (matches the publication title on line 30). |
| `docs/common_text/clair3_task.md` | Line 8 | Conditional comment '<!-- if: theiaviral-->' is missing the space before the closing '-->' that is used consistently elsewhere ('<!-- endif -->'). | For consistency, format as '<!-- if: theiaviral -->'. |
| `docs/common_text/gubbins_task.md` | Line 26 | Double space: "`tree_args`:  When" has two spaces after the colon. | Use a single space after the colon. |
| `docs/common_text/gubbins_task.md` | Line 27 | Double space: "alignment if  more than 25 %" has two spaces, and "25 %" has a space before the percent sign (inconsistent with "25%" style and unusual in American English). | Remove the extra space and write "25%". |
| `docs/common_text/ivar_trim_task.md` | Line 18 | Trailing whitespace after "aligned and sorted BAM file. " (line ends with a space). Also line 19 has trailing whitespace. | Remove trailing whitespace. |
| `docs/common_text/kaptive_task.md` | Lines 8 and 12 | Misspelling in the linked anchor text/URL context "baunannii" (should be "baumannii") appears in the Kaptive Wiki link anchors. While URLs are exempt, this reflects an upstream typo; the prose species name is correct. | No action needed if the URL anchor is correct upstream; otherwise note the species is _A. baumannii_, not "baunannii." |
| `docs/common_text/kaptive_task.md` | Lines 6, 7, 9, 11, 13 | Multiple lines have trailing whitespace at the end of paragraphs. | Remove trailing whitespace. |
| `docs/common_text/meningotype_task.md` | Line 6 | "_porB_ sequencing typing" is grammatically incorrect; should be "sequence typing" (consistent with the standard term and with the title's use of "Serotyping"). | Change "_porB_ sequencing typing" to "_porB_ sequence typing". |
| `docs/common_text/meningotype_task.md` | Line 5 (title) and Line 6 | The title says "Serotyping" but the body describes serogrouping and several other typing schemes, not serotyping; terminology is inconsistent. | Reconcile the title with the body (e.g., "Typing" or "Serogrouping & Typing") so the same concept is named consistently. |
| `docs/common_text/minimap2_task.md` | Line 29 | Missing period at the end of the sentence: "...please see the [Minimap2 manpage](...)" has no terminal punctuation. | Add a period after the link. |
| `docs/common_text/nanoq_task.md` | Line 14 (techdetails table, Original Publication row) | Missing closing table-cell pipe `\|` at the end of the Original Publication(s) row, unlike the other rows in the table. | Add a trailing `\|` after the publication link to match the table formatting of the other rows. |
| `docs/common_text/nextclade_task.md` | Line 11 | "_L. rabies_" appears to be an incorrect species abbreviation; rabies virus is _Lyssavirus rabies_ but "_L. rabies_" is non-standard/possibly wrong, and the linked dataset is for rabies. This may be a factual/typo error. | Verify and correct the species name (e.g., "rabies virus" or the correct italicized taxonomic name). |
| `docs/common_text/ngmaster_task.md` | Line 10 | Inconsistent terminology: "multi-antigen sequencing typing" and "sequencing typing" should be "sequence typing" (the standard term, also used as "sequence typing" in the title). | Change "sequencing typing" to "sequence typing" in both occurrences. |
| `docs/common_text/pangolin_task.md` | Line 8 | Misspelling: "Reoccuring" should be "Recurring" (note this is part of the scorpio acronym expansion). | Correct to "Recurring" (Serious Constellations of Reoccurring/Recurring Phylogenetically-Independent Origin). Verify the official scorpio expansion; the accepted spelling is "Reoccurring" with two r's at minimum, not "Reoccuring". |
| `docs/common_text/pirate_task.md` | Lines 6, 8 (trailing whitespace) | Trailing whitespace at end of several lines (minor formatting). | Remove trailing whitespace. Low priority. |
| `docs/common_text/plasmidfinder_task.md` | Line 19 | Trailing "<br>" with no following content after the Original Publication link. | Remove the dangling trailing "<br>" in the table cell. |
| `docs/guides/phylogenetics.md` | Line 84 (SNP threshold bullet list, lines 82-85) | The bulleted list under 'It can be difficult to determine SNP thresholds because of:' is indented with leading spaces under a paragraph, which in Markdown will likely render as a code block or fail to render as a proper list (the colon line is not a list item and the bullets are over-indented). | Remove the extra indentation so the dashes start at the left margin (or are nested correctly) to render as a real bulleted list. |
| `docs/workflows/genomic_characterization/freyja.md` | Line 34 (Freyja_Dashboard_PHB bullet) | Stray trailing '=' character after the link: '[**Freyja_Dashboard_PHB**](freyja.md#freyja_dashboard)='. | Remove the trailing '=' so the line ends with the link. |
| `docs/workflows/genomic_characterization/freyja.md` | Line 158 (freyja_demixed format table) | Mismatched quote in the example value: "[('Delta', 0.65), ('Other', 0.25), ('Alpha', 0.1')]" has a stray apostrophe after 0.1. | Change '0.1'' to '0.1' to fix the mismatched quote. |
| `docs/workflows/genomic_characterization/theiaviral.md` | Line 126 (heading) and Line 127 (`extract_unclassifed` ... set to `false`) | Within the `extract_unclassified` toggle block, the parameter name is spelled `extract_unclassifed` (missing the second 'i') in the body text on line 127, which is inconsistent with the heading and other references. | Correct the body text to `extract_unclassified` to match the parameter name used elsewhere. |
| `docs/workflows/phylogenetic_construction/ksnp3.md` | Line 16 (video link) | Inconsistent capitalization of the tool name: link text reads 'Using KSNP3' (all caps KSNP) while the rest of the page and the actual tool name use 'kSNP3'. | Change 'KSNP3' to 'kSNP3' in the link text for consistency. |
| `docs/workflows/standalone/kraken2.md` | Suggested databases table, 'viral w/human' row (line 41) | This row appears to be missing the Source column value: the cells are shifted so that '4.5' (size) lands in the Source column and the date/k-mer columns are off by one compared to every other row. The Source URL is missing. | Add the missing Source cell (or an empty cell) so the row's columns align with the header: Source, Database Size (GB), Date of Last Update, Bracken K-mer Lengths. |
| `docs/workflows/standalone/kraken2.md` | Line 93, Example Kraken2 report toggle | Spacing around the parenthetical em-dash phrase is inconsistent: 'very little -if any- read contamination' uses hyphens with surrounding spaces rather than proper dashes. | Use spaced en dashes or em dashes: 'very little — if any — read contamination'. |
| `docs/workflows/standalone/kraken2.md` | Lines 93, 96, 105 etc. | Inconsistent spacing before percent sign: '84.35 %' and '6 %' (space) versus '~2 %' — and elsewhere in the repo style is typically no space ('98.78%'). | Standardize to no space before the percent sign (e.g. '84.35%') to match the rest of the documentation. |
| `docs/workflows/standalone/ncbi_amrfinderplus.md` | Line 9 | Double space after the first sentence: 'stress genes.  Such AMR genes' and again 'point mutations.  In TheiaProk'. | Replace the double spaces with single spaces. |

## All Findings by File

### `docs/common_text/abricate_abaum_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | The term 'typing scheme' is used without explanation; a non-technical reader may not know what a typing scheme is or why a 95% identity match matters. | Briefly define 'typing scheme' (e.g., 'a standardized method for categorizing plasmids into types') so non-bioinformaticians understand the significance. |
| low | grammar | Line 6 | Inconsistent spacing in the identity threshold: ">/= 95 %" has a space before the percent sign, while the next paragraph writes "95%" without a space. | Standardize to ">= 95%" (no space before %), matching the "95%" used on line 8, and replace the ">/=" with the clearer ">=" or "at least 95%". |
| low | accessibility | Line 6 | Cross-reference link text 'the above section on Plasmid Identification' relies on the reader's position in a larger assembled page; in this standalone fragment the destination context is unclear. | Keep the descriptive link text but ensure the assembled page renders the anchor; consider 'the Plasmid Identification section above' to make the destination explicit. |

### `docs/common_text/abricate_bacterial_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 8 | The distinction between 'acquired resistance genes' and 'point mutations' is stated but not explained; a non-technical reader will not understand the practical consequence of ABRicate not detecting point mutations. | Add a short clarifying clause, e.g., 'It only detects acquired resistance genes (genes a bacterium has gained), NOT point mutations (small changes within existing genes), so some resistance may be missed.' |
| low | grammar | Line 8 | Missing serial (Oxford) comma in the database list: '...EcOH, PlasmidFinder, Ecoli_VF and VFDB.' | Add a comma before 'and': '...Ecoli_VF, and VFDB.' for consistency with American style used elsewhere. |
| low | ease-of-use | Line 8 | Acronyms AMR (implied), CARD, ARG-ANNOT, MEGARES, EcOH, VFDB appear without expansion; only VFDB is expanded later (line 10). | Since these are database proper names, expansion is optional, but consider noting that these are named databases of resistance/virulence genes for readers unfamiliar with them. |

### `docs/common_text/abricate_flu_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 8 | The sentence explaining that ABRicate 'typically works by screening contigs for the presence of acquired resistance genes' introduces a contrast that may confuse a non-technical reader, since 'contigs' is unexplained jargon. | Define 'contigs' on first use (e.g., 'contigs, the assembled stretches of DNA sequence') to aid comprehension. |

### `docs/common_text/abricate_vibrio_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 21 | Subject-verb agreement error: 'The presence or absence of specific genes are used...' The singular subject 'presence or absence' should take 'is'. | Change 'are used' to 'is used'. |
| low | accessibility | Line 8 (toggle heading 'Resistance Gene Database') | The toggle is labeled 'Resistance Gene Database', but the table contents are species markers, toxin genes, and serogroup markers, not resistance genes; the label is misleading. | Rename the toggle to something accurate like 'Vibrio Characterization Gene Database' so the label matches the table content. |
| low | ease-of-use | Lines 14-15 | Terms 'biotype' and 'serogroup' (and 'Classical'/'El Tor' biotypes, 'O1'/'O139' serogroups) are used without definition, which non-technical readers may not understand. | Add a brief parenthetical defining biotype and serogroup as ways of subclassifying the bacterium below the species level. |

### `docs/common_text/agrvate_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 6 | Redundant wording: 'Many S. aureus strains often have nonfunctional agr activity' uses both 'Many' and 'often', which is redundant. | Remove one, e.g., 'Many S. aureus strains have nonfunctional agr activity' or 'S. aureus strains often have...'. |
| low | ease-of-use | Line 6 | Terms 'locus', 'operon', and 'virulence phenotypes' appear without definition for a non-technical audience. | Add brief plain-language glosses, e.g., 'the agr locus (a specific region of the genome)' and 'virulence phenotypes (traits affecting how harmful the bacterium is)'. |
| low | ease-of-use | Line 8 | Jargon 'k-mers', 'in-silico PCR', and 'non-canonical agrD' is dense and unexplained, making the methods paragraph hard to follow for non-bioinformaticians. | Briefly define k-mers (short DNA sequence fragments) and in-silico PCR (a computational simulation of the lab PCR technique), or note that these are technical method details. |

### `docs/common_text/amr_search_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 16 | The phrase 'the system accounts for interactions between multiple SNPs, genes, and suppressors' uses the unexpanded acronym SNPs and the term 'suppressors' without definition. | Expand SNPs on first use ('single nucleotide polymorphisms (SNPs), single-letter changes in DNA') and briefly define 'suppressors' in this context. |
| low | accessibility | Line 16 | 'S/I/R classification' is defined inline, which is good, but the abbreviation appears bolded before the expansion making it momentarily unclear. | Acceptable as-is, but consider reordering to 'Sensitive/Intermediate/Resistant (S/I/R) classification' so the expansion precedes the abbreviation. |
| low | ease-of-use | Line 14 | 'in silico' is used without explanation; non-technical readers may not know this means 'by computer simulation'. | Add a brief gloss on first use, e.g., 'in silico (computational) antimicrobial resistance profiling'. |

### `docs/common_text/amrfinderplus_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 8 | This is a single very long, dense paragraph mixing what the tool does, taxa-specific behavior, intrinsic-gene filtering, and the GFF/protein FASTA option. It is hard to scan for a non-technical reader. | Break into shorter paragraphs or a bulleted list separating (1) what it detects, (2) taxa-specific results, and (3) the optional GFF/FASTA input. |
| low | ease-of-use | Line 8 | 'intrinsic genes', 'GFF', 'protein FASTA', and 'point mutations' are used without plain-language explanation. | Add brief glosses (e.g., 'intrinsic genes (genes naturally present in nearly all members of the species)') and define the GFF/FASTA file types or note they are genome annotation files. |
| low | ease-of-use | Line 10 | Long paragraph packs three reference resources and a database-update note together, which is dense. | Consider converting the three lookup resources into a short bulleted list for easier scanning. |

### `docs/common_text/arln_stats_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 8 | Metrics 'assembly ratio', 'percent GC statistics', and 'percent Q30' are listed without explanation of what they mean or why they matter for PASS/FAIL assessment. | Add brief definitions: assembly ratio (observed vs. expected genome size), GC content (proportion of G and C bases), and Q30 (the percentage of bases sequenced with high accuracy). |
| low | accessibility | Line 13 (toggle 'NCBI Assembly Stats Explained') | The toggle promises the columns are 'Explained' but only lists the column names verbatim without any explanation of what each represents. | Either rename the toggle to 'NCBI Assembly Stats Columns' or add a short description for each column (e.g., what StDev, CDS, and Consensus_TAXID mean). |
| low | ease-of-use | Line 10 | This is one dense paragraph describing the Q30 script, the assembly-stats file source, and the assembly-ratio calculation; it assumes familiarity with NCBI prokaryotic assembly data. | Split into separate short paragraphs (one per statistic source) and briefly note why the assembly ratio is compared against CDC reference data. |

### `docs/common_text/artic_consensus_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 8 | The stages description ('Alignment post-processing occurs, where primers are removed and various trimming steps are undertaken. Variants are detected...') uses passive, jargon-heavy phrasing that a non-technical reader will find hard to follow. | Reword in active voice and define key terms, e.g., 'reads are aligned to a reference genome; primer sequences are removed; then variants (differences from the reference) are identified and a consensus genome is built.' |
| low | ease-of-use | Line 10 | The Clair3 model default value is given but there is no guidance on how a non-technical user would know whether it is 'suitable for your sequencing data'. | Add a sentence on when to change it (e.g., 'choose a model matching your Nanopore basecalling configuration; consult your sequencing platform documentation if unsure'). |
| low | ease-of-use | Line 13 | Terms 'primer BED', 'BED formatting', and 'primer pools' appear without explanation in the info box. | Briefly note that a BED file is a plain-text file describing primer positions, for readers unfamiliar with the format. |

### `docs/common_text/artic_guppyplex_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 6 | 'chimeric reads' and 'amplicon-based viral sequencing' are unexplained jargon for a non-technical audience. | Add brief glosses, e.g., 'chimeric reads (artificial sequences formed by joining unrelated fragments)' and explain that amplicon-based sequencing targets specific regions using primers. |

### `docs/common_text/assembly_metrics_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 6 | 'BAM file', 'coverage', 'depth', 'average base quality', and 'average mapping quality' are listed without explanation; non-technical readers may not know these terms. | Add brief definitions for the key metrics (e.g., coverage = how much of the genome is represented; depth = how many reads cover each position) and note that a BAM file stores aligned sequencing reads. |

### `docs/common_text/augur_align_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 6 | The single-sentence description uses 'strips insertions relative to a reference sequence' without explaining why this is done, and assumes the reader knows what a multiple sequence alignment is. | Add a brief 'why' (e.g., 'aligning sequences makes them directly comparable for building phylogenetic trees') and gloss 'strips insertions' as removing positions not present in the reference. |

### `docs/common_text/augur_ancestral_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | Dense, jargon-heavy sentence with multiple undefined concepts (ancestral nucleotide sequences, maximum-likelihood, TreeTime, joint vs marginal models) that a non-technical reader cannot parse. | Add a plain-language sentence explaining what 'ancestral sequence reconstruction' achieves (e.g., 'estimates what the genetic sequences of common ancestors likely looked like') and briefly note that TreeTime is the underlying software. |
| low | accessibility | Line 6 / Line 8 (acronym) | Acronym 'TreeTime' and 'ML' concept appear without expansion of 'maximum-likelihood' context; also 'keep-ambiguous'/'infer-ambiguous' options on line 8 are referenced without explaining what 'ambiguous' bases are. | Briefly define ambiguous bases (e.g., positions where the nucleotide could not be determined) so the mutually-exclusive option note is meaningful to non-specialists. |
| low | ease-of-use | Line 6 | 'though "marginal" input is permitted' is ambiguous — it is unclear what 'marginal input' is or how a user would choose it. | Clarify what the joint vs marginal models do and which input controls the choice, or note this is an advanced option most users can leave at default. |

### `docs/common_text/augur_clades_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| high | typo | Line 6 | Malformed/escaped header string: 'header "clade\tgene\tsite\alt"' — the last separator is '\alt' (a single backslash) instead of '\t' (tab), inconsistent with the others. | Correct to '\talt' so the column separators are consistent (clade\tgene\tsite\talt). |
| medium | ease-of-use | Line 6 | Very dense single sentence describing the required clades_tsv format using technical terms ('tab-delimited', 'integer site', 'delineate clade-defining mutations') without a concrete example for a non-technical reader. | Break into shorter sentences and add a small example table or sample row showing the expected clades_tsv format. |
| low | ease-of-use | Line 6 | Term 'clades' and 'sequence signatures' are used in the heading and body without a plain-language definition of what a clade is. | Add a brief definition (e.g., 'a clade is a group of related samples that share a common ancestor and defining mutations'). |

### `docs/common_text/augur_export_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | grammar | Line 6 | 'can be inputted directly into' — 'inputted' is nonstandard; preferred form is 'input' or 'loaded into'. | Reword to 'can be loaded directly into the web-based tree viewer auspice.us' or 'can be uploaded directly to'. |
| low | ease-of-use | Line 6 | 'Auspice-formatted JSON' and 'Auspice' are introduced without explaining that Auspice is Nextstrain's interactive tree-visualization tool; a non-technical reader will not know what Auspice is. | Add a brief gloss, e.g., 'an Auspice-formatted JSON (Auspice is Nextstrain's interactive phylogenetic tree viewer)'. |

### `docs/common_text/augur_mutation_context_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | grammar | Line 6 | Subject-verb agreement error: 'These mutations have been shown to be a characteristic ... which indicate adaptation' — 'which' refers to the editing/mutations but uses 'indicate' where 'characteristic' is singular; the clause is also awkward. | Reword, e.g., 'These mutations are characteristic of APOBEC3-type editing, which indicates the virus is adapting to circulation among humans.' |
| medium | ease-of-use | Line 6 | 'APOBEC3-type editing' is highly specialized terminology presented without definition; non-technical readers will not understand what APOBEC3 is or why these mutations matter. | Add a brief plain-language explanation, e.g., 'APOBEC3 is a human enzyme that edits viral genomes; this mutation pattern is a signal that the virus has been adapting to spreading among humans.' |
| low | accessibility | Line 13 | 'NGA/TCN context of G→A or C→T mutations' uses nucleotide-context shorthand (N, NGA, TCN) without explaining that N denotes any nucleotide; unclear to non-specialists. | Add a short note that 'N' means any nucleotide and explain what the sequence context represents. |
| low | accessibility | Line 15 | Link text 'this example Nextstrain Mpox tree' is acceptable, but 'this' as a leading word is mildly vague; could be more descriptive. | Phrase the link text to fully describe the destination, e.g., 'an example Nextstrain Mpox tree showing these Color By options'. |

### `docs/common_text/augur_refine_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | Concepts 'calibrate branch lengths with respect to time', 'maximum-likelihood', and 'TimeTree' are unexplained for non-technical readers; also 'Augur_Prep_PHB' is referenced without context of where it comes from. | Add a plain-language note on what time-calibration produces (a tree where branch lengths reflect elapsed time) and clarify that Augur_Prep_PHB is an upstream task in the pipeline. |
| low | typo | Line 6 | Likely inconsistency: the surrounding Augur tasks describe time-calibration via 'TreeTime', but this file says 'TimeTree'. If the underlying software is TreeTime, this is a naming error. | Verify and correct the software name to 'TreeTime' if that is the tool used (the techdetails table on lines 12-14 references TreeTime, not TimeTree, confirming the body text is likely wrong). |
| low | grammar | Line 6 | Inconsistent capitalization: 'Time-Calibration' is capitalized mid-sentence (line 6) and also in the heading; mid-sentence it should be lowercase. | Use 'time-calibration' (lowercase) in the body sentence. |

### `docs/common_text/augur_traits_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | 'reconstruct the ancestral traits' is jargon; non-technical readers won't know what 'ancestral trait reconstruction' means or why they would want it. | Add a plain-language explanation, e.g., 'estimates likely characteristics (such as geographic location or lineage) of the ancestors in the tree.' |
| low | grammar | Line 6 | Slightly awkward phrasing 'to determine what trait metadata to use' — 'what' should be 'which' for a defined set of columns. | Reword to 'to specify which metadata columns to use as traits'. |

### `docs/common_text/augur_translate_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 6 | 'annotated features in the required reference_genbank input' assumes the reader knows what a GenBank file and 'annotated features' are. | Briefly note that the reference_genbank file describes where genes are located in the reference genome, which the tool uses to translate nucleotides into proteins. |

### `docs/common_text/augur_tree_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | Tradeoff between FastTree and IQ-TREE 2 ('performance improvements at the cost of phylogenetic resolution') and 'GTR evolutionary model' are stated without explaining the terms; a non-technical reader cannot judge when to switch tools or models. | Add brief plain-language guidance on when to choose FastTree vs IQ-TREE 2 and a one-line note on what a substitution/evolutionary model is. |
| low | grammar | Line 6 | Informal/awkward number style: 'if the alignment contains 100s of sequences' — '100s' is informal. | Use 'hundreds of sequences'. |

### `docs/common_text/bakta_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 18 | Missing article: 'stored in Google Cloud Storage Bucket' should be 'stored in a Google Cloud Storage bucket'. | Add the article 'a' and lowercase 'bucket': 'stored in a Google Cloud Storage bucket'. |
| low | accessibility | Line 5 (heading) | Heading says 'Assembly Annotation (alternative)' — 'alternative' is unexplained out of context (alternative to what?). | Clarify that Bakta is an alternative to the default annotation tool, e.g., note the default tool's name so 'alternative' has meaning. |
| low | ease-of-use | Line 8 | 'annotation' / 'regions of interest' are used without explaining what genome annotation produces for a non-technical reader. | Add a brief note that annotation identifies and labels genes and other features within the assembled genome. |

### `docs/common_text/bandage_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | Terms 'de novo assembly graphs', 'contigs', and 'misassemblies' are unexplained; the non-technical reader will not know what a contig or assembly graph is. | Add brief definitions, e.g., 'contigs are continuous stretches of assembled sequence' and 'an assembly graph shows how these pieces connect'. |

### `docs/common_text/bbduk_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 14 | Double space: 'with [BBDuk](...)  (*"Bestus...' has two spaces before the parenthetical. | Remove the extra space so there is a single space before the parenthetical. |
| low | grammar | Line 19 | Missing comma after introductory word: 'Additionally primers will be trimmed' should be 'Additionally, primers will be trimmed'. | Add a comma after 'Additionally'. |
| low | accessibility | Line 12 | 'PhiX' is introduced in the heading and line 12; while it is later explained ('PhiX is a viral genome...'), the heading uses the acronym/term before definition. | This is minor since the body defines it; consider ensuring the first body mention defines PhiX before referencing it as a 31-mer match target. |
| low | ease-of-use | Lines 12, 19, 23 | Technical terms '31-mer', 'kmer size', and 'degenerate primers' are used without definition for a non-technical audience. | Add a brief gloss for k-mer (a short subsequence of length k used for matching) on first use, and a note on what degenerate primers are. |

### `docs/common_text/bbmap_reformat_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | Terms 'interleaved FASTQ files', 'deinterleaves', and 'repairs them' are unexplained; non-technical readers won't know what interleaving/deinterleaving means. | Add a brief explanation that interleaved files combine paired reads into one file, and deinterleaving separates them back into two files (read1 and read2). |
| low | grammar | Line 9 | Trailing whitespace and slightly awkward phrasing: 'given the necessity to return paired end reads.' Also 'paired end' should be hyphenated as 'paired-end'. | Hyphenate 'paired-end' and reword, e.g., 'because it must return paired-end reads.' |
| low | ease-of-use | Line 9 | 'IRMA' acronym is used without expansion or definition for the TheiaCov audience. | Expand IRMA on first use (Iterative Refinement Meta-Assembler) or link to its documentation. |

### `docs/common_text/bcftools_consensus_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 8 | Acronyms/jargon 'VCF', 'indels', and 'FASTA' are used without expansion; the non-technical target audience may not know these terms. | On first use, expand: 'VCF (Variant Call Format) file', 'indels (insertions and deletions)', and 'FASTA format (a standard text format for storing sequences)'. |
| medium | ease-of-use | Line 8 | This is a single dense sentence packing the whole workflow (filtering, aligning, normalizing, indexing, generating, substituting, masking). Hard for a non-technical reader to follow. | Break the description into a short intro sentence plus a bulleted list of the steps the task performs. |
| low | grammar | Line 8 | Phrase 'the `min_depth` and `min_allele_freq` input parameter' uses singular 'parameter' for two parameters. | Change to 'input parameters' (plural) to agree with the two parameters listed. |
| low | grammar | Line 8 | 'left aligns' is missing a hyphen when used as a verb phrase paralleling other verbs; also reads slightly awkwardly in the run-on list. | Use 'left-aligns and normalizes indels' for consistency, and consider breaking the long sentence into two for readability. |

### `docs/common_text/busco_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 7 | British spelling 'Orthologue' in the acronym expansion (and 'orthologs' American is used in the publication title on line 28). Inconsistent spelling within the file. | The acronym does spell BUSCO with 'Orthologue', so this is constrained; if changing, note that BUSCO officially uses 'Orthologs'. Consider 'Orthologs' for American-English consistency, but verify against official BUSCO naming first. |
| low | grammar | Line 20 | Missing serial comma: 'low duplicated (D), fragmented (F) and missing (M) percentages' lacks the Oxford comma used elsewhere in the doc set. | Add comma: 'duplicated (D), fragmented (F), and missing (M) percentages'. |
| low | ease-of-use | Line 11 | The term 'BUSCO notation' string is shown before the abbreviations are defined; a non-technical reader sees a cryptic string first. | Add a brief lead-in such as 'BUSCO summarizes results in a compact notation. For example:' before showing the string. |

### `docs/common_text/bwa_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | typo | Lines 8, 11, 15, 19, 22 (all 'Burrow-Wheeler Aligner') | BWA is consistently expanded as 'Burrow-Wheeler Aligner'; the correct name is 'Burrows-Wheeler Aligner' (with an 's'). | Correct all occurrences to 'Burrows-Wheeler Aligner' (matches the publication title on line 30). |
| low | ease-of-use | Line 8 | 'assembly_fasta' is a raw variable/file name dropped into prose without explanation for a non-technical reader. | Phrase as 'the Pilon-polished assembly (assembly_fasta)' or describe it as the assembled genome file. |
| low | ease-of-use | Line 22 (theiaviral) | 'skani' and 'samtools' appear without context for what they are. | Briefly note these are tools, e.g. 'a reference genome either selected by skani (a sequence comparison tool) or provided by the user'. |
| low | ease-of-use | Lines 8, 11, 15, 19, 22 | 'BAM file' is used repeatedly without expansion; non-technical readers may not know it is an alignment file. | On first use expand: 'BAM file (a compressed file storing read alignments)'. |

### `docs/common_text/cauris_cladetyper.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 24 | Trailing whitespace after the sentence; also 'a close relative annotation' reads awkwardly. | Remove trailing space and reword to 'an annotation from a close relative, such as Clade IV'. |
| low | accessibility | Lines 14-21 (clade reference table) | The table is presented with only a lead-in line ('See more information ... below') and no textual explanation of what the columns mean (e.g., what Genome Accession or BioSample Accession represent) for non-technical readers or screen-reader users. | Add a sentence describing the table's purpose and briefly defining the accession columns, e.g. 'Each clade is represented by a publicly available reference genome, identified by its NCBI Genome and BioSample accession numbers.' |
| low | ease-of-use | Line 7 | 'GAMBIT' is used as an acronym without expansion on first use (the full expansion appears only in the publication title on line 32). | Expand on first use: 'GAMBIT (Genomic Approximation Method for Bacterial Identification and Tracking)'. |
| low | ease-of-use | Line 10 | 'genomic signature comparison' is jargon unlikely to be understood by the non-technical audience. | Add a brief plain-language gloss, e.g. 'compared against this database by matching characteristic patterns in the genome (genomic signatures)'. |

### `docs/common_text/centroid_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 8 | 'mash distances' is unexplained jargon; a non-technical reader will not know what mash is or what a mash distance measures. | Add a brief gloss, e.g. 'pairwise Mash distances (a fast estimate of how genetically similar two genomes are)'. |
| low | grammar | Line 8 | Double space after 'generate the tree.' and trailing whitespace at end of line. | Use a single space between sentences and remove trailing whitespace. |
| low | ease-of-use | Line 8 | 'most central genome' / 'most central' is intuitive-but-undefined; the doc does not explain why centrality matters. | Briefly explain centrality, e.g. 'the genome that is, on average, most similar to all the others (the centroid)'. |

### `docs/common_text/cg_pipeline_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 5 (heading) | Heading 'Assessment of Read Quality, and Estimation of Genome Coverage' has an unnecessary comma before 'and' joining two noun phrases. | Remove the comma: 'Assessment of Read Quality and Estimation of Genome Coverage'. |
| low | ease-of-use | Line 6 | 'QUAST' acronym/tool name appears without expansion or explanation. | Briefly identify QUAST, e.g. 'the estimated genome length generated by QUAST (a genome assembly evaluation tool)'. |
| low | ease-of-use | Line 6 | 'coverage' (genome coverage) is a technical concept used without definition for non-technical readers. | Add a one-line definition of coverage, e.g. 'coverage (the average number of times each position in the genome was sequenced)'. |

### `docs/common_text/checkv_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | Jargon 'integrated proviruses', 'genome fragments', and 'closed genomes' are used without definition; difficult for non-technical public health scientists. | Add brief glosses, e.g. 'proviruses (viral genomes integrated into a host genome)' and 'closed genomes (complete, gap-free genomes)'. |
| low | grammar | Line 6 | Parallelism: 'identification of host contamination ..., estimating completeness ..., and identification of closed genomes' mixes noun and gerund forms. | Make parallel: 'identifying host contamination ..., estimating completeness ..., and identifying closed genomes'. |
| low | ease-of-use | Line 8 | 'weighted_contamination' and 'weighted_completeness' output fields are described tersely; the concept of 'weighted by contig length' may be unclear. | Add a short explanation, e.g. 'these averages give more influence to longer contigs (continuous stretches of assembled sequence).' |

### `docs/common_text/clair3_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | 'pileup-based calling' and 'full-alignment analysis' are technical terms used without explanation for the non-technical audience. | Briefly gloss these, or simplify to plain language describing a two-pass approach for finding genetic differences. |
| low | typo | Line 8 | Conditional comment '<!-- if: theiaviral-->' is missing the space before the closing '-->' that is used consistently elsewhere ('<!-- endif -->'). | For consistency, format as '<!-- if: theiaviral -->'. |
| low | grammar | Line 30 | Inconsistent capitalization/spelling of 'vcf' vs 'VCF': 'the VCF file' (line 30) but 'A filtered vcf file' (line 30). | Use 'VCF' consistently. |
| low | ease-of-use | Line 9 | 'ONT' acronym is used without expansion. | Expand on first use: 'ONT (Oxford Nanopore Technologies) data'. |
| low | ease-of-use | Line 16 | 'basecaller' is jargon used without definition. | Add a brief gloss, e.g. 'depending on the basecaller (the software that converts raw sequencer signal into DNA letters) and data type'. |

### `docs/common_text/clockwork_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 14 (techdetails table) | Publication title 'cohorts of bacterial genome' should be 'genomes' (plural) — verify against the actual Minos paper title. | Confirm and correct to 'cohorts of bacterial genomes' if the original title is plural. |
| low | ease-of-use | Line 5 (heading) | 'Illumina PE' uses the abbreviation 'PE' without expansion; non-technical readers may not know it means paired-end. | Expand to 'Illumina paired-end (PE) only' for clarity. |
| low | ease-of-use | Line 6 | 'H37Rv reference genome', 'BWA', and 'TBProfiler' are used without context; 'TB' is also assumed. | Briefly note H37Rv is the standard M. tuberculosis reference genome, identify BWA as the alignment tool, and define 'TB' as tuberculosis on first use. |

### `docs/common_text/concatenate_illumina_lanes_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | accessibility | Line 13 (techdetails table) | Table row is missing the trailing pipe '\|' present in other techdetails tables; minor formatting inconsistency that can affect rendering/screen readers. | Add the closing '\|' after the link for consistent table formatting. |
| low | ease-of-use | Line 6 | Assumes the reader understands multi-lane Illumina sequencing and the read1/read2 lane naming; no 'why' context for when this applies. | Add a sentence explaining when this is relevant, e.g. 'When a sample is sequenced across multiple lanes on an Illumina instrument, the FASTQ files must be combined before analysis.' |

### `docs/common_text/concatenate_variants_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | accessibility | Lines 20-24 (variant table) | The example table uses many undefined column abbreviations (CHROM, POS, TYPE, REF, ALT, EVIDENCE, FTYPE, NT_POS, AA_POS, LOCUS_TAG) with no textual key, making it inaccessible to non-technical readers and screen-reader users. | Add a short legend defining the key columns (e.g., CHROM = chromosome/contig, POS = position, REF/ALT = reference/alternate base, etc.). |
| low | ease-of-use | Line 16 | References 'the cat_files task' which a reader of this fragment may not have encountered; comparison may not aid understanding. | Briefly describe what cat_files does, or rephrase to avoid relying on prior knowledge of that task. |

### `docs/common_text/consensus_qc_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 6 | 'degenerate bases' is jargon used without definition for the non-technical audience. | Add a brief gloss, e.g. 'degenerate bases (ambiguous positions where more than one base is possible)'. |
| low | ease-of-use | Line 6 | '"N" bases' is used without explanation of what an N represents. | Note that 'N' bases are unknown or unresolved positions in the sequence. |

### `docs/common_text/contaminant_check_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | grammar | Line 6 | "based on if contaminant/host sequences pass thresholds" uses informal "based on if" construction. | Rewrite as "based on whether contaminant/host sequences pass thresholds". |
| medium | ease-of-use | Line 6-8 | Acronyms/terms 'breadth of coverage', 'depth of coverage', and 'FASTA' are used without explanation for a non-technical audience; thresholds min_expected_seq/max_unexpected_seq are referenced without saying what they control. | Briefly define breadth vs depth of coverage and note that FASTA is a sequence file format; add a one-line description of what min_expected_seq and max_unexpected_seq set. |
| medium | ease-of-use | Line 10 | "reported in JSON mappings" assumes the reader knows what a JSON mapping is. | Add a brief gloss, e.g. "reported as JSON (a structured text format) listing each sequence with its statistics." |
| low | grammar | Line 6 | Repeated awkward use of "inputted" / "previously inputted" reads stiffly and "previously inputted" is redundant. | Use "provided" or "input" e.g. "a comma-delimited string of expected_sequences" and "Each sequence from the contaminant/host FASTA". |
| low | grammar | Line 8 | "is reported depicting which expected_sequences failed" — "depicting" is the wrong word for text describing reasons. | Use "is reported indicating which expected_sequences failed and why". |
| low | ease-of-use | Whole task | No 'why/when to use this' guidance — the reader is not told the purpose (detecting contamination or host carryover before downstream analysis). | Add an opening sentence explaining when and why a user would run a contaminant check. |

### `docs/common_text/cophylogeny_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | Terms 'phylogeny', 'tips', 'branch', and 'topological' may be unfamiliar; "topological (tip/branch arrangement)" partially defines it but assumes phylogeny tree vocabulary. | Add a brief note that a phylogeny is an evolutionary tree, and that 'tips' are the samples at the ends of branches. |
| low | grammar | Line 6 | "two inputted phylogenies" — awkward use of "inputted". | Use "two input phylogenies" or "two provided phylogenies". |
| low | ease-of-use | Line 8 | Dense single paragraph covering two output files plus a caveat about evolutionary distance; hard to scan. | Split into a short bulleted list of the two outputs (with/without branch lengths) and a separate note about interpreting evolutionary distance. |

### `docs/common_text/data_summary_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| high | ease-of-use | Line 12 | Very long, dense paragraph mixing input format, behavior, output, and an optional setting; difficult for a non-technical reader to follow. | Break into bullets: (1) what to enter in data_summary_* and sample_names, (2) input format example, (3) what the task produces (CSV), (4) the phandango_coloring option. |
| medium | ease-of-use | Line 12 | 'CSV', 'Phandango', 'Newick tree', and 'AMR' (in the figure caption) are used without definition for a non-technical audience. | Expand on first use: CSV (comma-separated values file), Phandango (a tree/metadata visualization tool), Newick (a tree file format), AMR (antimicrobial resistance). |
| low | accessibility | Line 25 | Link text is the bare URL <http://jameshadfield.github.io/phandango/#/main> rather than descriptive text. | Use descriptive link text such as "the Phandango web viewer" pointing to the URL. |
| low | ease-of-use | Line 12 | "indicates presence (TRUE) or absence (empty)" but the example CSV on lines 17-21 shows FALSE values, which is inconsistent with the description. | Reconcile the description with the example (explain when FALSE vs empty appears), or correct one to match the other. |

### `docs/common_text/digger_denovo_wf.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | The example given is about fungal genomes ('Assembly of fungal genomes from short-reads...'), but the next paragraph says digger_denovo is used in TheiaProk and TheiaEuk; the fungal-only example may confuse readers about scope. | Generalize the example to apply to the relevant organisms (bacterial/eukaryotic), or clarify why a fungal example is used here. |
| medium | ease-of-use | Line 6 | Terms 'contigs', 'contiguous sequence', and 'short-reads' are used without definition; 'de novo' is italicized but never plainly defined for non-technical readers. | Add brief glosses: contig (a continuous stretch of assembled sequence), short reads (the fragments produced by Illumina sequencing). |
| low | grammar | Line 6 | "the process or product of attempting to reconstruct" is wordy and slightly ambiguous. | Simplify to "the process of reconstructing a genome from scratch ... using sequence reads." |
| low | ease-of-use | Line 8 | 'subworkflow', 'Shovill', and 'SPAdes' are referenced; the in-joke about the name 'digger' adds no functional value for the target audience. | Optional: keep the joke but ensure Shovill/SPAdes are clearly labeled as assembly software. |

### `docs/common_text/dnaapler_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | 'contigs', 'reorients', 'CDS', and gene names (dnaA, terL, repA, COG1474) are used with no explanation; CDS in particular (used on lines 14-16) is never expanded. | Expand CDS on first use (coding sequence) and add a one-line explanation of why reorienting a contig's start point matters. |
| low | grammar | Lines 9-16 | Inconsistent terminal punctuation in the bullet list — some bullets end with a period (lines 9, 14) and others do not (lines 10, 11, 12, 13, 15, 16). | Make terminal punctuation consistent across all bullets. |
| low | ease-of-use | Whole task | No 'why/when to use this' framing — the reader isn't told why assembly orientation matters or which mode to pick for their organism. | Add a short intro sentence on the purpose of reorientation and a pointer that mode should match the molecule type (chromosome/plasmid/phage/archaea). |

### `docs/common_text/download_terra_table_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 6 | 'TSV' and 'Terra project/workspace/table' are used without explanation; TSV is not expanded for a non-technical reader. | Expand TSV (tab-separated values file) on first use. |
| low | ease-of-use | Line 6 | No 'why/when to use this' guidance for a one-sentence task. | Add a brief note on when a user would want to download a Terra table (e.g., to archive or use it outside Terra). |

### `docs/common_text/ectyper_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | 'serotyping', 'pathotype', 'Shiga toxin typing', and 'in silico' are specialized terms; only pathotype is partially explained. | Add brief glosses for serotyping (classifying strains by surface antigens) and note 'in silico' means computational/from sequence data. |
| low | ease-of-use | Line 6 | Dense paragraph combining purpose, capabilities, and a definition of pathotype. | Split the pathotype definition into its own sentence or bullet for readability. |

### `docs/common_text/emmtyper_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 8 | 'BLAST', 'CDS', 'genome assembly', 'alleles', 'in silico PCR', and 'emm cluster' are used without definition for a non-technical audience. | Briefly expand/define BLAST (a sequence-comparison tool) and 'in silico PCR' (a simulated PCR from sequence data); gloss allele. |
| low | ease-of-use | Line 6 vs 8 | Tool name capitalization is inconsistent: title and code use 'emmtyper'; the prose also uses 'emmtyper' but the related file uses 'emm-typing-tool' — fine here, but ensure 'emmtyper' casing matches consistently (it does). Minor: 'M protein' vs 'M-type' vs 'M-type' usage could confuse. | Optionally clarify the relationship between the M protein, emm gene, and the reported emm-type. |

### `docs/common_text/emmtypingtool_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | accessibility | Line 8 | Link text "decision tree" points to a .png image; screen-reader users following the link land on an image with no described content, and the figure itself has no textual explanation in the doc. | Note that the link opens an image, and/or summarize what the decision algorithm does in text. |
| medium | ease-of-use | Line 8 | 'maps the reads', 'bowtie2', 'alleles', '% coverage', '% identity' used without definition for a non-technical audience. | Briefly gloss read mapping and percent identity/coverage; identify bowtie2 as a read-alignment tool. |
| low | ease-of-use | Line 5 | The '==for Illumina_PE only==' highlight uses 'Illumina_PE' which non-technical readers may not parse (paired-end Illumina). | Expand to 'Illumina paired-end (PE) reads only'. |

### `docs/common_text/estimate_genome_length_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | 'taxon', 'accession', 'completeness status', and 'basepairs' are used without definition; 'top reference accession' is unclear to a non-technical reader. | Briefly define taxon (a named group of organisms), accession (a database identifier), and basepairs; clarify what 'top reference accession' means. |
| low | grammar | Line 6 | "an inputted taxon" — awkward use of "inputted". | Use "a given taxon" or "the input taxon". |
| low | ease-of-use | Line 6 | Long single sentence chaining multiple actions (acquire metadata, retrieve accession, generate summary, calculate average). | Split into shorter sentences or a short bulleted list of what the task produces. |

### `docs/common_text/ete4_identify_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | 'taxonomy hierarchy', 'taxonomic rank', and the downstream module names are technical; non-technical readers may not know what a 'taxon ID' or 'taxonomic rank' is. | Add a brief explanation that taxonomic rank means levels like species, genus, family, etc. (which is partly shown later on line 12). |
| low | grammar | Line 6 | "a user's inputted taxonomy" — awkward; also slightly redundant with 'user's'. | Use "a user-provided taxon" or "the input taxon". |
| low | grammar | Line 11 | Missing period after 'a.k.a' in the toggle title "`rank` a.k.a `read_extraction_rank`". | Use 'a.k.a.' with the trailing period, or 'also known as'. |

### `docs/common_text/extract_clade_mutations_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| high | ease-of-use | Line 6 | Highly technical with no plain-language framing: 'Augur-compatible clades.tsv', 'signature clade-defining sequences', 'nucleotide JSON outputted by Augur Ancestral', 'amino acid JSON outputted by Augur Translate' will be opaque to non-technical readers. | Add an opening sentence in plain language explaining the purpose (identifying the mutations that define each clade/group in a tree) before the technical details, and gloss Augur, JSON, clade. |
| medium | ease-of-use | Line 8 | 'monophyletic clades', 'mutation signatures', and 'clade metadata column' assume phylogenetics knowledge. | Briefly define monophyletic clade (a group containing an ancestor and all its descendants) and clarify what the clade metadata column is and where it comes from. |
| low | grammar | Line 6 | "A nucleotide JSON outputted by Augur Ancestral is required" — 'outputted' reads awkwardly (twice on lines 6). | Use "produced by" or "output by" consistently. |

### `docs/common_text/fasta_utilities_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 6 | "fasta" is written in lowercase here but elsewhere file formats are capitalized (e.g., FASTA/FASTQ). Inconsistent casing of the format name. | Use "FASTA" to match the conventional capitalization used elsewhere in the docs. |
| low | ease-of-use | Line 6 | "index a reference fasta file" assumes the reader knows what indexing a FASTA file means and why it is necessary; a non-technical user gets no context for the purpose of this step. | Add a short clause explaining why, e.g. "...uses samtools to index a reference FASTA file, which creates a companion file that lets downstream tools quickly look up positions in the reference." |

### `docs/common_text/fastp_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 9 (read_qc_trim block) | This is a single very long, dense sentence packed with multiple parameters, defaults, and parenthetical exceptions (window size 4, 10 for TheiaEuk, trim_quality_trim_score defaults, trim_minlen defaults). It is hard for a non-technical reader to parse. | Break into shorter sentences or a small bulleted list of parameters and their defaults (window size, quality score, minimum length) for paired-end vs single-end. |
| low | grammar | Line 9 | Unbalanced/nested parentheses: "a sliding window (with a default window size of 4 (10 for TheiaEuk), specified with `trim_window_size`)" stacks parentheses inside parentheses, making the sentence hard to follow. | Rephrase to avoid nested parentheses, e.g. "...a sliding window whose default size is 4 (or 10 for TheiaEuk), set with `trim_window_size`." |
| low | ease-of-use | Line 9 and Line 24 | "sliding window" is used as a technical term without explanation; a non-bioinformatician may not know what trimming with a sliding window means. | Add a brief plain-language gloss the first time, e.g. "a sliding window (a small section that moves across the read)". |
| low | ease-of-use | Line 16 (Default read-trimming parameters table) | Table rows -5 20 and -3 20 both have the identical explanation "enables read end-trimming", which is confusing because the flags differ (5' vs 3' end) but the explanations do not distinguish them. | Differentiate the explanations, e.g. "enables trimming from the 5' (front) end" and "enables trimming from the 3' (tail) end". |
| low | ease-of-use | Line 11 | References "Trimmomatic's default configuration" without context; a reader who has not used Trimmomatic will not know what comparison is being drawn or why it matters. | Briefly note that Trimmomatic is the alternative/default trimming tool being compared against, so the comparison is meaningful. |

### `docs/common_text/fastplong_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | Long, dense sentence combining window size, quality score, minimum length, and directionality (cut_front/cut_tail) into one block. Hard for a non-technical reader to follow. | Split into shorter sentences or bullets; separate the trimming defaults from the directionality explanation. |
| low | grammar | Line 8 | Awkward phrasing "a FASTA, a string of start, or a string of end adapters can be specified" — "a string of start" reads as if a word is missing. | Rephrase, e.g. "a FASTA file of adapters, a start-adapter string, or an end-adapter string can be specified with `fastplong_adapter_fasta`, `fastplong_start_adapter`, or `fastplong_end_adapter` respectively." |
| low | accessibility | Line 5 (heading) and body | "ONT" appears in the heading "ONT Read Trimming" and is never expanded anywhere in this fragment; first-use acronym is undefined for non-technical readers. | Expand on first use, e.g. "Oxford Nanopore Technologies (ONT) read trimming". |

### `docs/common_text/fastq_scan_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 9 | "If QC has been performed correctly" uses the acronym QC without expansion; non-technical readers may not know QC = quality control. | Expand on first use: "If quality control (QC) has been performed correctly". |

### `docs/common_text/fastqc_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 5 (heading) and Line 8 | The heading calls FastQC "Read Quantification (alternative)" and line 8 describes it as quantifying reads (text copied from fastq-scan), but FastQC's main purpose is read quality assessment/visualization, mentioned only as an afterthought on line 10. This mislabels the tool's primary function and could confuse users choosing between tools. | Reframe so the primary purpose (quality assessment and visualization) is stated first, and clarify how/why a user would choose FastQC over fastq-scan. |
| low | accessibility | Line 16 | Inconsistent capitalization of "GitHub" — link text reads "FastQC on Github" (lowercase h) while every other file in this set uses "GitHub". | Capitalize as "FastQC on GitHub" for consistency. |
| low | ease-of-use | Line 8 | "If QC has been performed correctly" uses the acronym QC without expansion. | Expand on first use: "If quality control (QC) has been performed correctly". |

### `docs/common_text/filter_contigs_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | accessibility | Lines 16-18 (digger block) | The bulleted list of skip parameters is indented with extra spaces inside an already-indented admonition; rendered Markdown may not display these as a proper list. The surrounding list items use 4-space indentation but these use 8. | Verify the rendered output and align list indentation so the parameter bullets render as a list rather than a code/preformatted block. |
| low | ease-of-use | Line 7 / Line 12 | "homopolymer contigs (contigs of any length that consist of a single nucleotide)" uses the jargon "homopolymer" and "contigs"; while contigs is glossed by context, homopolymer is only partially explained and may still be unclear to a non-technical reader. | Slightly expand, e.g. "homopolymer contigs (sequences made up of a single repeated nucleotide, such as AAAAA), which are usually sequencing artifacts." |
| low | ease-of-use | Line 23 | "contigs that meet specified criteria" is vague — it does not restate which criteria (length, coverage, homopolymer) so the sentence adds little for a reader skimming. | Reference the criteria explicitly, e.g. "...only contigs that meet the length, coverage, and homopolymer criteria described above." |

### `docs/common_text/flu_antiviral_substitutions_wf.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | grammar | Line 6 | Broken/ambiguous sentence: "...of H1N1 or H3N2 flu samples, or any in non-subtype-specific PA, PB1, and PB2 segments." The phrase "or any in" is grammatically incomplete and unclear. | Rephrase, e.g. "...in the HA, NA, and MP segments of H1N1 or H3N2 flu samples, or in the non-subtype-specific PA, PB1, and PB2 segments." |
| medium | accessibility | Line 6 | Multiple flu segment acronyms (HA, NA, MP, PA, PB1, PB2) and subtypes (H1N1, H3N2) are used with no expansion or explanation; a non-technical reader will not know these are influenza gene segments. | Add a brief gloss on first use, e.g. "the HA (hemagglutinin), NA (neuraminidase), and MP (matrix protein) segments" and note PA/PB1/PB2 are polymerase segments. |
| low | accessibility | Line 8 | Link text "a list of known amino-acid substitutions associated with antiviral resistance" points to a .wdl task file on GitHub (source code), which is unexpected; the descriptive text implies a readable list but the destination is code. Slightly misleading link target. | Either clarify the link goes to the source code containing the list, or point to a more human-readable rendering of the list. |
| low | ease-of-use | Line 12 | The format example "`Protein:AAPositionAA`" is terse; a non-technical reader may not map "AAPositionAA" to "original amino acid, position, new amino acid". | Spell out the format, e.g. "in the format `Protein:[original amino acid][position][new amino acid]`, for example `NA:V95A`". |

### `docs/common_text/flye_denovo_wf.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 9 (tip box) | The Medaka model guidance is dense and assumes familiarity with terms like "basecaller model", "pore", "caller variant", and the model-string format. A non-technical reader is unlikely to know how to determine their basecaller model. | Add brief guidance on where a user can find their basecaller/model information (e.g., from the sequencing run report or with their sequencing core) and define "basecaller" in plain terms. |
| low | grammar | Line 11 | Missing comma after introductory clause: "If Flye is being run on legacy data the Medaka model will likely be...". | Add a comma: "If Flye is being run on legacy data, the Medaka model will likely be...". |
| low | grammar | Line 13 | Awkward passive phrasing "Recently generated data will likely be suited by the default model". | Rephrase to active voice, e.g. "The default model `r1041_e82_400bps_sup_v5.0.0` is likely suitable for recently generated data." |
| low | accessibility | Line 6 and Line 9 | "ONT" is used without expansion in this fragment; first-use acronym undefined for non-technical readers. | Expand on first use: "Oxford Nanopore Technologies (ONT) data". |

### `docs/common_text/flye_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | "repeat graphs", "de Bruijn graphs", and "exact k-mer matches" are advanced bioinformatics concepts presented without definition; non-technical readers will not understand the comparison. | Either simplify to the practical takeaway (Flye tolerates the higher error rate of long reads) or add brief plain-language definitions of these terms. |
| low | accessibility | Line 6 and elsewhere | "ONT" is used repeatedly without expansion in this fragment. | Expand on first use: "Oxford Nanopore Technologies (ONT)". |
| low | ease-of-use | Lines 17-22 (flye_read_type table) | Terms such as "Guppy5+ SUP", "Q20", "CLR", and "HiFi" appear without expansion; a non-technical user selecting a read type may not recognize these basecaller/chemistry terms. | Add brief expansions on first use (e.g., CLR = Continuous Long Read, HiFi = High-Fidelity) or a short note pointing to where users can identify their read type. |
| low | ease-of-use | Line 9 (theiaviral block) | Slightly confusing: "It can be enabled by setting the `call_raven` parameter to `false`" — enabling Flye by setting a parameter named after a different tool (Raven) to false is non-obvious and may confuse readers. | Clarify the relationship, e.g. "Flye is used instead of Raven when `call_raven` is set to `false`." |

### `docs/common_text/gambit_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | accessibility | Line 6 | "k-mer based approach" uses the term k-mer without definition; non-technical readers will not know what a k-mer is. | Add a brief gloss, e.g. "a k-mer based approach (comparing short DNA sub-sequences of fixed length)". |

### `docs/common_text/gamma_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 8 | Dense paragraph comparing GAMMA to AMRFinder and BLAT vs BLAST without explaining what AMRFinder, BLAT, or BLAST are; non-technical readers cannot evaluate the comparison. | Briefly note BLAT and BLAST are sequence-alignment search tools, and that AMRFinder is an alternative AMR-detection tool, so the comparison is meaningful. |
| low | grammar | Line 8 | Trailing double space and run-on style; also "a protein identity based tool" should be hyphenated as a compound modifier. | Use "a protein-identity-based tool" and remove trailing whitespace. |
| low | accessibility | Line 5 (heading) and body | "AMR" appears in the heading "AMR Genotyping (optional)" and is never expanded; first-use acronym undefined. | Expand on first use: "antimicrobial resistance (AMR) genotyping". |
| low | ease-of-use | Line 10 | "multifasta database", "coding sequences", "hypervirulence", "plasmid markers", and "prokaryotic database" are stacked jargon terms; combined with "GFF output" and "--gff flag", the paragraph assumes substantial domain knowledge. | Define or simplify the key terms (especially "coding sequences" and "GFF") or link to a glossary for non-technical readers. |

### `docs/common_text/gene_coverage_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | accessibility | Lines 30, 32-37, 41-45 (BED file sections) | "BED file" / "BEDfiles" and "GBFF" are used extensively without expansion or explanation of what these file formats are; non-technical readers will not know these formats or how to create one. | Expand acronyms on first use (BED = a tab-delimited file of genomic regions; GBFF = GenBank Flat File) and briefly describe what each contains. |
| low | grammar | Line 30 | Inconsistent spelling of the file format: "BEDfiles" (one word) here versus "BED file" (two words) on lines 32-37 and 42-45. | Standardize to "BED files" throughout. |
| low | accessibility | Lines 8-24 (JSON map examples) | The two JSON code-block examples are shown with no textual explanation of what the numbers mean (e.g., that 99.5 is a percentage and 10 is a depth value), so a screen-reader user or non-technical reader cannot interpret them. | Add a sentence noting that values in percent_coverage_by_gene are percentages and values in depth_by_gene are average read depths. |
| low | ease-of-use | Line 5 (heading) and Line 6 | "Depth and Breadth of Coverage" are domain terms; while line 6 mentions "average read depth" and "percent of a region covered", the terms depth and breadth are not explicitly tied to those definitions for a non-technical reader. | Briefly define: depth = how many reads cover each position on average; breadth = the percentage of the region covered above the minimum depth. |

### `docs/common_text/genoflu_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | Jargon not defined for non-technical readers: "whole-genome genotype," "segment," and "clade 2.3.4.4b" are used without explanation. A non-bioinformatician may not know flu genomes are divided into segments or what a clade is. | Add a brief plain-language note, e.g. that influenza genomes are made up of 8 separate RNA segments, and that a "clade" is a related group of viruses sharing a common ancestor. |
| low | grammar | Line 6 | Wrong article: "a H5N1" should use "an" because "H" is pronounced "aitch" (vowel sound). Appears twice ("a H5N1 (...) flu sample"). | Change "a H5N1" to "an H5N1" in both occurrences. |

### `docs/common_text/genotyphi_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | grammar | Line 6 | Software name appears inconsistent/likely a typo: "SeqSero2S" is referenced in prose, but the variable is "seqsero2s_predicted_serotype". Elsewhere this tool is commonly "SeqSero2". Also the phrase "identification of the \"Typhi\" serotype by SISTR or SeqSero2S" pairs a tool list awkwardly. | Verify the correct tool name (SeqSero2 vs SeqSero2S) and use it consistently with the variable name. |
| medium | ease-of-use | Line 6-8 | Acronyms/jargon undefined for non-technical readers: SNV (line 8), "quinolone-resistance determining regions," "plasmid replicons," "k-mer based genotyping," "lineages, clades, and subclades." | Expand SNV (single nucleotide variant) on first use and add brief definitions for k-mer and clade/lineage, or link to a glossary. |
| low | grammar | Line 6 | Subject-verb / relative clause: "subtypes of the IncHI1 plasmid which is associated with multidrug resistance" needs a comma before "which" for the non-restrictive clause. | Change to "the IncHI1 plasmid, which is associated with multidrug resistance." |

### `docs/common_text/gubbins_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | grammar | Line 15 | Inconsistent capitalization and naming of the tool: "gubbins" (lowercase) on line 15 vs "Gubbins" elsewhere. Also "RaxML" (lines 18, 26) is misspelled; the correct tool name is "RAxML". | Capitalize "Gubbins" consistently and correct "RaxML" to "RAxML". |
| medium | ease-of-use | Line 18, Line 20 | Undefined jargon for non-technical readers: SNP/SNPs, "loci," "phylogenies," "relaxed selection," "mutational hotspots," "recombinant regions." Also "industry standard" is informal. | Expand SNP (single nucleotide polymorphism) on first use and add brief plain-language definitions or a glossary link for loci/phylogeny. |
| low | typo | Line 26 | Double space: "`tree_args`:  When" has two spaces after the colon. | Use a single space after the colon. |
| low | typo | Line 27 | Double space: "alignment if  more than 25 %" has two spaces, and "25 %" has a space before the percent sign (inconsistent with "25%" style and unusual in American English). | Remove the extra space and write "25%". |
| low | grammar | Line 18 and Line 26 | Inconsistent model naming: "GTR-GAMMA" (line 18), "GTR+GAMMA" (line 15), and "GTRGAMMA" (line 26) refer to the same substitution model with three different spellings. | Standardize the substitution model name across the file. |

### `docs/common_text/hicap_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6, Line 8 | Undefined jargon for non-technical readers: "cap locus," "serotype/serology," "in silico," "post-translation modification." The audience may not know what a locus is or what serotyping means. | Add a one-sentence plain-language intro before the quote explaining serotyping (classifying strains by surface antigens) and that a locus is a region/location of genes on the genome. |
| low | grammar | Line 8 (quoted block) | British spelling "characterisation" should be American "characterization". (Note: this is inside a quoted block from upstream docs; flag only if house style requires American throughout.) | If house style mandates American spelling, change to "characterization" (or note it is a direct quote). |

### `docs/common_text/iqtree1_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6, Line 8 | Undefined jargon for non-technical readers: "maximum-likelihood method," "nucleotide substitution model," "bootstrapping," "clade," "SH-aLRT." No explanation of why these support thresholds (80%/95%) matter. | Add a brief plain-language note on what a phylogeny/clade is and that bootstrap/SH-aLRT values measure how confident we can be in a branch of the tree. |

### `docs/common_text/iqtree2_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6, Line 12 | Undefined jargon for non-technical readers: "Model Finder," "maximum-likelihood," "nucleotide substitution model," "bootstrapping," "clade," "SH-aLRT," "core_genome." | Add brief definitions or a glossary link; explain what the support-value thresholds mean in plain language. |
| low | ease-of-use | Line 8 | Ambiguous wording: "whether any sites were masked with recombination" reads as if sites are masked _with_ recombination; the intent is sites masked _because of_ recombination. | Reword to "whether any sites were masked due to recombination." |

### `docs/common_text/irma_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | Acronym LABEL is used without expansion or explanation, and "typing and subtyping," "segments," and "meta-assembler" are unexplained for non-technical readers. | Briefly explain that LABEL is the read-sorting component, and define typing/subtyping in plain language. |
| low | grammar | Line 6 | Awkward repetition and phrasing: "please review the ABRicate task outputs which will provide this information" repeats "this information" twice in adjacent sentences ("For determining this information... which will provide this information"). | Rephrase, e.g. "To distinguish Victoria from Yamagata lineages, review the ABRicate task outputs." |
| low | grammar | Line 6 | Comma splice / run-on: "iteratively maps reads to a collection of reference sequences ... and iteratively edits the references to account for high population diversity and mutational rates that are characteristic of Influenza genomes" is a very long sentence; also "mutational rates" reads better as "mutation rates." | Split the sentence and consider "mutation rates." |
| low | ease-of-use | Line 6 | The "Note:" about Victoria/Yamagata is embedded at the end of a dense paragraph and easy to miss. | Pull the Note out into a separate admonition or its own line for visibility. |

### `docs/common_text/ivar_consensus_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | Undefined jargon for non-technical readers: "reference-based consensus assembly," "allele frequency," "depth," "stringency." No "why/when to use this" guidance. | Add a one-line plain-language description of what a consensus assembly is and what these thresholds control. |
| low | ease-of-use | Lines 20, 23, 26 | Inconsistent terminology: parameter prose says "minimum mapping quality" / "minimum read depth" / "minimum allele frequency" but the intro (line 6) says "minimum quality," "minimum depth," "minimum allele frequency." Same concepts, different wording across the same fragment. | Use consistent terms (e.g. always "minimum mapping quality" and "minimum read depth"). |
| low | ease-of-use | Lines 22-23, 25-26 | The min_map_quality and min_allele_freq parameter descriptions both say they affect "variant calling and subsequent consensus sequence generation," but this is the consensus task; mapping-quality/allele-freq wording is copied from the variants task and may confuse readers about what this task does. | Verify the descriptions match the consensus task's actual behavior and adjust wording so it is not identical to the variants fragment. |

### `docs/common_text/ivar_trim_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | grammar | Line 18 | Wrong article and broken phrasing: "a [_organism-specific parameters_-determined]" should be "an" before a vowel sound; the embedded link label makes the sentence hard to parse. | Reword, e.g. "Using the user-provided (or organism-specific) `primer_bed` file..." with the link applied to a cleaner phrase, and use "an" if a vowel sound follows. |
| medium | ease-of-use | Line 18, Line 22 | Undefined jargon/acronyms for non-technical readers: "soft-clips," "BAM file," "BED format," "0-based index," "sliding window approach," "base quality." These are dense for the target audience. | Add brief definitions (e.g. BAM = aligned reads file; BED = a simple text format listing coordinates) or link to a glossary. |
| low | typo | Line 18 | Trailing whitespace after "aligned and sorted BAM file. " (line ends with a space). Also line 19 has trailing whitespace. | Remove trailing whitespace. |
| low | ease-of-use | Line 22 | Ambiguous parenthetical: "a sliding window approach (default: 4)" does not state what the 4 refers to (window size in bases). | Clarify, e.g. "a sliding window approach (default window size: 4 bases)." |

### `docs/common_text/ivar_variants_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6, Line 8 | Undefined jargon for non-technical readers: "samtools mpileup," SNVs/indels (expanded but still technical), VCF, "allele frequencies," "stringency," "intermediate variants." | Briefly define VCF (a standard text file listing genetic variants) and add a plain-language note on why intermediate-frequency variants matter. |
| low | ease-of-use | Lines 22, 25 | Parameter descriptions for min_map_quality and min_allele_freq say they affect "consensus sequence generation," but this is the variant-calling task; the consensus reference may confuse readers about the task's purpose. | Verify wording matches the variants task and adjust so it describes variant calling specifically. |

### `docs/common_text/kaptive_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | grammar | Line 6 | Subject-verb agreement error: "Specificity for non-antibiotic therapeutics (e.g. phage therapy) bear particular significance" — the subject "Specificity" is singular and requires "bears." | Change "bear" to "bears." |
| medium | ease-of-use | Line 6, Line 8, Line 10 | Very dense paragraphs with heavy undefined jargon for non-technical readers: "capsular polysaccharide (CPS)," "epidemiological marker," "K locus (KL)," "O antigen," "lipooligosaccharide," "outer core (OC) locus," "serotype," "virulence." Several long sentences pack multiple concepts. | Break the long paragraphs into shorter ones and add brief plain-language definitions for capsule/serotype/locus, or a glossary link, so non-specialists can follow why this typing matters. |
| low | typo | Lines 8 and 12 | Misspelling in the linked anchor text/URL context "baunannii" (should be "baumannii") appears in the Kaptive Wiki link anchors. While URLs are exempt, this reflects an upstream typo; the prose species name is correct. | No action needed if the URL anchor is correct upstream; otherwise note the species is _A. baumannii_, not "baunannii." |
| low | typo | Lines 6, 7, 9, 11, 13 | Multiple lines have trailing whitespace at the end of paragraphs. | Remove trailing whitespace. |
| low | grammar | Line 6 | "e.g." used without a following comma (American style typically uses "e.g.,"). Occurs on line 6 ("e.g. phage therapy"). | Add a comma: "e.g., phage therapy." |

### `docs/common_text/kleborate_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 5 (task title) / Line 6 | The acronyms MLST and AMR appear in the title and body with no expansion. A non-technical reader will not know MLST means Multi-Locus Sequence Typing or that AMR means antimicrobial resistance. | Expand on first use, e.g. 'MLST (Multi-Locus Sequence Typing)' and 'AMR (antimicrobial resistance)'. |
| low | ease-of-use | Line 6 | Terms 'serotype', 'capsular (K antigen)', 'lipopolysaccharide (LPS) (O antigen)', and 'ICE_Kp_' are presented without plain-language context for a non-bioinformatician. | Add a brief clause explaining what serotyping/K and O antigens indicate, or link to a glossary. |
| low | ease-of-use | Line 8 | 'K antigen and O antigen locus typing via _wzi_ alleles' assumes specialist knowledge with no plain explanation of what Kaptive adds or why a user would want it. | Add a short 'when/why' note explaining the benefit of enabling Kaptive. |

### `docs/common_text/kmerfinder_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 8 | The core concept of '_k_-mers' (16-mers, 11-mers) is used heavily without a one-line definition for a non-technical reader. | Add a brief definition, e.g. 'k-mers (short sub-sequences of length k from the DNA)'. |
| low | ease-of-use | Line 8 | 'prokaryotic species' and 'coding regions' are unexpanded jargon for a low-technical audience. | Briefly gloss 'prokaryotic (e.g. bacteria)' on first use. |
| low | ease-of-use | Line 8 | 'training data' is referenced but the reference database was the subject of the sentence; the switch in terminology (database vs training data) may confuse readers about whether these are the same thing. | Use consistent terminology, e.g. 'reference database' throughout instead of switching to 'training data'. |

### `docs/common_text/kraken2_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | grammar | Line 46 | Broken/non-parallel list item: 'inputted in place of Kraken reports in downstream tasks' lacks a subject/verb consistent with the other bullets, which start with verbs like 'increases', 'decreases'. It reads as a sentence fragment. | Rephrase for parallelism, e.g. 'is used in place of Kraken reports in downstream tasks, such as `qc_check` and `krona`'. |
| medium | ease-of-use | Line 40 | 'Bayesian model to probabilistically estimate read abundances' is dense jargon for a non-technical audience with no plain explanation of what Bracken practically does. | Add a plain-language summary, e.g. 'Bracken statistically re-estimates how many reads belong to each species to give more reliable counts.' |
| low | grammar | Line 47 | 'outputted separate of the' is awkward/ungrammatical; 'separate of' should be 'separately from'. | Change to 'is output separately from the `kraken/kraken2_report`'. |
| low | grammar | Lines 40, 23, 32 | Inconsistent representation of the boolean: text says set to "true"/"false" with quotation marks in some places (lines 23, 32, 40) while other fragments use backticks (`false`). Inconsistent formatting for the same concept. | Standardize boolean values to backticks (`true`/`false`) across all fragments. |
| low | accessibility | Line 30, 48, 49 | Trailing whitespace / stray blank lines with spaces (lines 30, 48 end with spaces); minor but can affect rendering consistency. | Remove trailing whitespace. |
| low | ease-of-use | Line 49 | Very dense paragraph mixing k-mer database selection, naming conventions, and skip conditions; hard for a non-expert to follow. | Break into short bullets or sentences separating (1) default behavior, (2) how to override, (3) when Bracken is skipped. |

### `docs/common_text/kraken_parser_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | 'scatter portion of the workflow' and 'scatter shards' are pipeline-engineering jargon (WDL scatter) that a non-technical public health scientist will not understand, and they are never defined. | Briefly explain in plain terms, e.g. 'the part of the workflow that runs many parallel jobs (one per taxon)'. |
| low | ease-of-use | Line 6 | 'lightens the computation load' / 'lowers the number of scatter shards' describes internal optimization but offers no 'why this matters to me' for the user. | Add a short note on the user-facing benefit (e.g. faster, cheaper runs) or clarify this is an internal step requiring no user action. |

### `docs/common_text/krakentools_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | 'taxon ID' is used without explanation; a non-technical reader may not know this is an NCBI numeric identifier for an organism. | Briefly gloss, e.g. 'taxon ID (the NCBI numeric identifier for an organism)'. |
| low | ease-of-use | Line 13 | 'descendant taxa' and 'subordinate taxonomic levels' assume understanding of taxonomic hierarchy; could confuse the target audience. | Add a brief plain-language example (e.g. 'all organisms grouped beneath that taxon, such as all species within a genus'). |

### `docs/common_text/krona_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 6 | 'hierarchical datasets' and 'multi-layered pie charts' may be unclear; the sentence is fairly dense for a non-technical reader. | Add a brief plain-language framing, e.g. 'Krona creates interactive, zoomable pie charts that show which organisms were found and how they relate to each other.' |

### `docs/common_text/ksnp_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Lines 6/9 (title) and Line 11 | 'pan-genome and core-genome phylogenies' and 'Phylogenetic Construction' are unexpanded specialist terms; no plain explanation of what this workflow produces or when to use it. | Add a short sentence explaining the output (a tree showing how samples are related) and define pan-genome vs core-genome briefly. |
| low | ease-of-use | Line 11 | Missing 'why/when to use this' guidance for the workflow, which the audience needs. | Add one sentence on when a user would choose kSNP (e.g. building SNP-based phylogenetic trees from assemblies without a reference genome). |

### `docs/common_text/legsta_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 5 (title) / Line 6 | SBT (sequence-based typing) is the central concept but the acronym is only loosely tied to the spelled-out form; a non-technical reader may not connect 'SBT' on later uses. | On first use write 'sequence-based typing (SBT)' explicitly, then use SBT. |
| low | accessibility | Line 8 | Trailing whitespace after 'Technical Details' heading line (line 8 ends with spaces). | Remove trailing whitespace. |
| low | ease-of-use | Line 6 | '_in silico_' is used without explanation for a non-technical audience. | Gloss as '_in silico_ (computationally predicted)'. |

### `docs/common_text/lissero_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | accessibility | Lines 12-29 (Serogroup Decision Tree code block) | The decision tree is presented as a code block with no textual/plain-language explanation of how to read it; screen-reader users and non-programmers may struggle with the Python-like if/elif/else logic. | Add a short sentence explaining how to interpret the tree (gene presence/absence determines serogroup) or summarize the outcomes in prose/a table. |
| low | accessibility | Line 8 | Trailing whitespace at end of line 8. | Remove trailing whitespace. |
| low | ease-of-use | Line 6 | 'serogroup' and 'somatic (O) or flagellar (H) biosynthesis' are specialist terms used without plain-language explanation. | Briefly explain what a serogroup is, or link to a glossary for the antigen terminology. |

### `docs/common_text/mapping_stats_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 6 | 'BAM file', 'coverage', 'depth', 'average base quality', 'average mapping quality' are technical metrics presented without explanation for a non-technical audience. | Briefly define BAM (aligned-reads file) and add a one-line note on what coverage/depth indicate, or link to a glossary. |
| low | ease-of-use | Line 6 | 'reported on a per sequence basis' is slightly ambiguous (per reference sequence? per read?). | Clarify, e.g. 'reported separately for each reference sequence/contig'. |

### `docs/common_text/medaka_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 5 (title) and body | 'Polishing', 'Flye assembly', 'consensus sequence', and 'basecaller' are used without plain-language definitions for a non-technical reader. | Add a brief gloss for 'polishing' (correcting errors in the draft assembly) on first use. |
| low | grammar | Line 14 | 'If both automatic model selection fails and no user-specified model is provided' — 'both ... and' is awkward here since the two clauses are not both subjects of 'fails'. | Rephrase, e.g. 'If automatic model selection fails and no user-specified model is provided'. |
| low | ease-of-use | Lines 12-14 | The 'model' concept (line 10 'the model that was used to generate the read data') is central but never explained in plain terms; a non-technical user will not know what a Medaka model is or how to find theirs. | Add a brief sentence explaining that the model corresponds to the sequencing chemistry/basecaller settings, and that automatic selection usually handles this. |

### `docs/common_text/megahit_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 13/15 | '_De novo_ assembly', 'k-mer strategy', 'succinct de Bruijn graphs', and 'metagenomics' are specialist terms with no plain-language explanation for a non-technical audience. | Add a brief gloss for '_de novo_ assembly' (building a genome from reads without a reference) and avoid/explain 'de Bruijn graphs'. |
| low | ease-of-use | Line 18 | Logical inconsistency that may confuse users: the text says the task is 'turned off by default' and 'only called if MetaviralSPAdes fails', then says 'It can be enabled by setting the `skip_metaviralspades` parameter to `true`' — enabling a fallback by a 'skip' parameter is counterintuitive and the relationship between the two statements is unclear. | Clarify the trigger logic: explain whether MEGAHIT runs automatically on SPAdes failure or only when `skip_metaviralspades` is true, and reconcile the two sentences. |

### `docs/common_text/meningotype_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | typo | Line 6 | "_porB_ sequencing typing" is grammatically incorrect; should be "sequence typing" (consistent with the standard term and with the title's use of "Serotyping"). | Change "_porB_ sequencing typing" to "_porB_ sequence typing". |
| medium | ease-of-use | Line 6 | Multiple acronyms/terms are undefined for a non-technical audience: "_in silico_", "finetyping", "serogrouping", and "BAST" (BAST is partially expanded but the concept is dense). The whole capability list is a single dense sentence. | Add a one-line plain-language definition of "_in silico_" (computational/from sequence data) and briefly explain what serogrouping/finetyping accomplish, or break the function list into a bulleted list for readability. |
| low | typo | Line 5 (title) and Line 6 | The title says "Serotyping" but the body describes serogrouping and several other typing schemes, not serotyping; terminology is inconsistent. | Reconcile the title with the body (e.g., "Typing" or "Serogrouping & Typing") so the same concept is named consistently. |
| low | ease-of-use | Line 8 | "BLAST" and "PubMLST" are introduced without expansion or explanation, which a non-technical reader may not recognize. | Briefly note that BLAST is a sequence-comparison algorithm and PubMLST is a public reference database of allele sequences. |

### `docs/common_text/metabuli_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | Dense, jargon-heavy paragraph for the target audience: "k-mer structure", "metamer", "AA conservation", "sensitive homology detection", "DNA mutations for specific differentiation between closely related taxa" are unexplained. | Add a plain-language sentence on what read classification does (assigning each sequencing read to the organism it came from) and briefly gloss "k-mer" and "homology" for non-bioinformaticians. |
| low | grammar | Line 31 | "the memory and cpus allocated" mixes a full word with an abbreviation and uses lowercase "cpus" inconsistently with the parameter naming. | Use "the memory and CPUs allocated" or refer to the `memory`/`cpu` parameters consistently. |
| low | ease-of-use | Line 28 | "taxonkit-generated taxdump", "taxonomy hierarchy", and abbreviations GTDB/NCBI are used without explanation; a non-technical reader will not know what a taxdump file is. | Add a brief note that a taxdump file is the reference taxonomy (the naming/classification system) Metabuli uses to label reads. |

### `docs/common_text/microreact_export_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 6 | "API", "access token", "Terra tables", and "JSON" are used without expansion for a non-technical audience. | Expand "API" (application programming interface) on first use and briefly explain that an access token is a personal credential used to authenticate with Microreact. |
| low | ease-of-use | Line 8 | "a response JSON will be returned to the user for review" does not explain what the user should look for or do with it. | Briefly state what the response JSON indicates (e.g., whether the submission succeeded or what errors occurred) so the user knows how to act on it. |

### `docs/common_text/midas_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | accessibility | Lines 17-23 (MIDAS report table) | The example table is presented before its columns are explained; a screen-reader user encounters the raw data table with no preceding textual description of what it represents. | Add a brief sentence before the table summarizing what it shows (one row per detected species with read counts, coverage, and relative abundance), in addition to the column descriptions that follow. |
| medium | ease-of-use | Line 32 | Very long, dense paragraph explaining how four output columns are derived; hard to follow for the target audience. | Break the derivation of each output column (`midas_primary_genus`, `midas_secondary_genus`, `midas_secondary_genus_abundance`, `midas_secondary_genus_coverage`) into a bulleted list. |
| low | ease-of-use | Line 10 | Acronyms "WGS" and "QC" and the term "non-target taxa" / "relative frequency" appear without expansion for a non-technical reader. "co-opted" is informal jargon. | Expand "WGS" (whole-genome sequencing) and "QC" (quality control) on first use, and consider replacing "co-opted" with "adapted". |
| low | ease-of-use | Lines 40, 44, 45 | "ANI", "SNP", "pan-genome", "non-redundant genes", and "single-copy genes" are technical terms used without definition. | Expand "ANI" (average nucleotide identity) and "SNP" (single-nucleotide polymorphism) on first use and briefly gloss pan-genome. |
| low | ease-of-use | Line 57 | "Google Cloud Storage bucket" and "Terra workspace" instructions assume familiarity; non-technical users may not know how to obtain/provide the link. | Add a brief pointer or example of what the `midas_db` link should look like (e.g., a gs:// path) to reduce ambiguity. |

### `docs/common_text/minimap2_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | accessibility | Lines 5-9 (conditional title block and heading hierarchy) | The `??? task` admonition header is inside an `if` conditional, so when neither `long_read_flags` nor `only_map_ont` is set, the indented body content (line 9 onward) has no parent header. This can produce orphaned/unheaded content for some workflow builds. | Verify every rendering path provides a task header so body text is never orphaned under a missing heading. |
| low | typo | Line 29 | Missing period at the end of the sentence: "...please see the [Minimap2 manpage](...)" has no terminal punctuation. | Add a period after the link. |
| low | ease-of-use | Lines 9, 14, 18, 22, 26 | Acronyms/formats "ONT", "SAM", "PAF", and "aligner" are used without expansion, and raw parameter strings are shown without explaining the user need not understand them. | Expand "ONT" (Oxford Nanopore Technologies) on first use and briefly note that SAM/PAF are standard alignment output file formats; add a note that the listed flags are preset defaults the user does not need to set. |
| low | ease-of-use | Line 9 | "'modes' are a group of preset options" is vague for a non-technical reader and does not explain why mode choice matters to them. | Add a sentence clarifying that the workflow selects the appropriate mode automatically based on the data type, so no user action is required. |

### `docs/common_text/mummer_ani_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 10 | "between 2 genomes" uses a numeral for a small number mid-sentence, inconsistent with surrounding prose style. | Write out "two genomes". |
| low | ease-of-use | Line 8 | "taxa" usage: "the sequence is the same taxa as the reference" uses the plural "taxa" with a singular subject. | Change to "the same taxon as the reference". |
| low | ease-of-use | Lines 10, 12 | Acronyms "CDC EDLB", "PulseNet", and "perl" capitalization, plus "enteric pathogens", may be unfamiliar to a non-technical reader. | Expand "CDC EDLB" (CDC Enteric Diseases Laboratory Branch) on first use and briefly note PulseNet is a national surveillance network; capitalize "Perl". |
| low | ease-of-use | Lines 14-17 | The two thresholds and their default values are introduced without explaining the units (percent) or what a user should do if a sample falls below them. | Clarify the values are percentages and add brief guidance on interpreting results that do not surpass both thresholds. |

### `docs/common_text/nanoplot_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 6 | "Read Quantification" title vs. body describing quality scores and read lengths; also "QC" is unexpanded. The 'why fewer clean reads than raw' is good but "if QC has been performed correctly" assumes the reader knows what QC is. | Expand "QC" (quality control) on first use; consider noting that read filtering/trimming removes low-quality reads, which is why fewer remain. |
| low | ease-of-use | Line 9 | Dense single sentence with nested parenthetical explaining why NanoPlot runs outside read_QC_trim_ont; hard to parse for a non-technical reader. | Split into two sentences and plainly state the purpose: it computes read statistics so accurate coverage can be estimated using the actual assembly length when no genome length is provided. |

### `docs/common_text/nanoq_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | typo | Line 14 (techdetails table, Original Publication row) | Missing closing table-cell pipe `\|` at the end of the Original Publication(s) row, unlike the other rows in the table. | Add a trailing `\|` after the publication link to match the table formatting of the other rows. |
| low | ease-of-use | Line 6 | "quality scores lower than 10" does not explain the scale (Phred quality), which a non-technical reader will not understand. | Briefly note these are Phred quality scores and what a higher value means, or link to a definition. |

### `docs/common_text/ncbi_datasets_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 10 | Long sentence with nested parenthetical ("(for all other genome types)") makes the virus-vs-genome package logic hard to follow. | Rephrase to clearly state: the virus package is used for viral assemblies and the genome package for all other types. |
| low | ease-of-use | Lines 16-17 | "gbff file" and "fasta file" formats are referenced without explanation; non-technical readers may not know the difference or implications. | Briefly note that a gbff file includes gene annotations while a fasta file is sequence only, which is why annotations may differ. |

### `docs/common_text/ncbi_scrub_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 8 | Highly technical single sentence: "k-mer database constructed of k-mers from Eukaryota derived from all human RefSeq records with any k-mers found in non-Eukaryota RefSeq records subtracted from the database" is very dense and uses undefined terms (k-mer, Eukaryota, RefSeq). | Break into shorter sentences and add a plain-language summary: the tool recognizes human DNA fragments and removes reads matching them, after excluding fragments shared with non-human organisms. |
| low | ease-of-use | Line 5 (title) and Line 6 | "HRRT" is used as the title acronym before it is expanded; the expansion appears in line 6 but the title-first encounter is unexplained. | Consider expanding HRRT in or near the title, or ensure the body expansion is the first encounter. |

### `docs/common_text/nextclade_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| high | typo | Line 11 | "_L. rabies_" appears to be an incorrect species abbreviation; rabies virus is _Lyssavirus rabies_ but "_L. rabies_" is non-standard/possibly wrong, and the linked dataset is for rabies. This may be a factual/typo error. | Verify and correct the species name (e.g., "rabies virus" or the correct italicized taxonomic name). |
| low | ease-of-use | Lines 6, 8 | Terms "clade", "mutation calling", "phylogenetic placement", "node", and "reference tree" are used without definition for a non-technical audience. | Add a brief plain-language definition of "clade" (a group of related viral genomes sharing a common ancestor) on first use. |
| low | ease-of-use | Line 6 | "calling variants between the two sequences" uses "variants" / "calling" jargon without explanation. | Briefly gloss that this means identifying the differences (mutations) between the sample and the reference. |

### `docs/common_text/ngmaster_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Lines 6-8 | Several acronyms are used without expansion: "NG-MAST", "WGS", "PCR", "NG-STAR", "23S rRNA". "NG-MAST" and "NG-STAR" are only loosely expanded later (line 10), after first use. | Expand WGS (whole-genome sequencing) and PCR (polymerase chain reaction) on first use, and move the NG-MAST / NG-STAR expansions to their first appearance. |
| low | typo | Line 10 | Inconsistent terminology: "multi-antigen sequencing typing" and "sequencing typing" should be "sequence typing" (the standard term, also used as "sequence typing" in the title). | Change "sequencing typing" to "sequence typing" in both occurrences. |
| low | ease-of-use | Line 6 | The paragraph order is confusing: NG-MAST is described first, then NG-STAR (line 8), then line 10 reveals NGMASTER combines both. A reader may be confused about the relationship until the end. | Lead with the line 10 summary (NGMASTER combines NG-MAST and NG-STAR) so the reader has context before the detailed descriptions. |

### `docs/common_text/organism_parameters_wf.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 (intro sentence) | "acquires and propagates default files and variable values for specific organisms" is jargon-heavy and abstract for a non-technical reader. It is unclear what 'propagates' means or why the user would care. | Rephrase in plainer terms, e.g. "This workflow automatically selects and applies the correct reference files and default settings for the organism you are analyzing, so you do not have to set them by hand." |
| low | accessibility | Line 5 and toggles (MPXV, RSV-A, RSV-B) | Acronyms MPXV, RSV, HA, NA, VOC are used without expansion on first use. A non-technical reader may not know MPXV = monkeypox/mpox virus, RSV = respiratory syncytial virus, or that HA/NA are influenza gene segments. | Expand acronyms on first use, e.g. "MPXV (mpox virus)", "RSV-A (respiratory syncytial virus A)", and note that HA and NA refer to the hemagglutinin and neuraminidase gene segments. |
| low | ease-of-use | Lines 28-60 (Flu toggle, H1N1/H3N2 etc.) | The flu defaults list reference files under bare labels 'HA' and 'NA' without explaining these are influenza gene segments, which non-technical users may not recognize. | Add a brief note that HA and NA denote the hemagglutinin and neuraminidase gene segments used for typing. |

### `docs/common_text/pangolin_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | typo | Line 8 | Misspelling: "Reoccuring" should be "Recurring" (note this is part of the scorpio acronym expansion). | Correct to "Recurring" (Serious Constellations of Reoccurring/Recurring Phylogenetically-Independent Origin). Verify the official scorpio expansion; the accepted spelling is "Reoccurring" with two r's at minimum, not "Reoccuring". |
| medium | ease-of-use | Line 8 | Dense single paragraph packed with terms (hash, designation cache, constellations, pangoLEARN, UShER, inference pipeline) that a non-technical reader will struggle to follow. | Break into shorter sentences or a brief bulleted list of steps, and add a one-line plain-language summary of what the user receives (a lineage assignment for each sample). |
| low | accessibility | Line 8 | Acronyms VOC and QC appear without expansion (VOC is expanded as 'variant of concern' but QC is not). | Expand QC as "quality control (QC)" on first use. |

### `docs/common_text/pasty_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 8 | "The serogroup can be predicted based off of those results." — "based off of" is colloquial/nonstandard. | Change to "based on those results." |

### `docs/common_text/pbptyper_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 6 | "The Penicillin-binding proteins (PBP) are responsible..." — the definite article plus plural reads awkwardly, and PBP is given as singular for a plural noun. | Rephrase to "Penicillin-binding proteins (PBPs) are responsible..." and use PBPs for the plural. |
| low | grammar | Line 6 vs Line 5/8 | Inconsistent capitalization of the tool name: "PBPTyper" (line 6) vs "PBPtyper" (lines 5, 8, 10). | Use one consistent spelling of the tool name throughout; standardize to "PBPtyper". |
| low | accessibility | Line 6 | Acronym ANI is later expanded in line 8 but beta-lactam and MIC are introduced; MIC is expanded but the relationship is dense. Acceptable, minor. | Optional: keep MIC expansion; no change strictly required. |

### `docs/common_text/phylovalidate_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| high | ease-of-use | Line 6 and Line 8 | The description is highly technical (polytomies, topologies, Lin-Rajan-Moret distance, matching cluster distance, Robinson-Foulds distance) with no plain-language explanation of what the tool does for the user or when they would use it. The audience is explicitly non-technical. | Add an opening plain-language sentence describing the purpose, e.g. "PhyloValidate checks whether two phylogenetic trees (evolutionary relationship diagrams) are essentially the same shape, which is useful for confirming that results are reproducible." Then keep the technical detail for advanced users. |
| medium | ease-of-use | Line 6 | "validate if the distance between these two phylogenies' topologies is less than an inputted `max_distance` float (0 by default)" assumes the reader knows what a topology distance and a 'float' are. | Briefly define topology (the branching shape of the tree) and avoid the term 'float' for this audience, or note it means a decimal number. |
| low | grammar | Line 6 | Subject-verb agreement: "Trees can only be compared if the number of nodes between the trees are the same." — subject "number" is singular. | Change "are the same" to "is the same". |
| low | grammar | Line 6 | "the tips must be the same between trees" combined with the following clause is a long run-on sentence. | Split into two sentences for readability. |

### `docs/common_text/pilon_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | accessibility | Line 9, Line 14 | Acronym BWA (line 9) and tool name metaSPAdes are used without context; BWA is an aligner not defined. | Optionally note BWA is a read alignment tool on first mention. |
| low | ease-of-use | Line 9 (if: digger block) | "small indels" uses unexpanded jargon (indels = insertions/deletions). | Expand on first use: "small indels (insertions and deletions)". |

### `docs/common_text/pirate_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | accessibility | Line 6 | Multiple unexpanded acronyms/terms for a non-technical audience: GFF3, CDS, orthologous gene families, pangenome, polytomies-adjacent concepts. CDS in particular is unexpanded. | Expand CDS as "coding sequences (CDS)" and briefly define pangenome and orthologous gene families on first use. |
| medium | ease-of-use | Lines 6-8 | Both paragraphs are dense, comma-spliced lists describing parameters and outputs; hard for non-technical readers to parse what they get. | Convert the list of generated outputs (line 8) into a bulleted list and simplify the parameter aside in line 6. |
| low | typo | Lines 6, 8 (trailing whitespace) | Trailing whitespace at end of several lines (minor formatting). | Remove trailing whitespace. Low priority. |

### `docs/common_text/plasmidfinder_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | typo | Line 19 | Trailing "<br>" with no following content after the Original Publication link. | Remove the dangling trailing "<br>" in the table cell. |

### `docs/common_text/polypolish_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | accessibility | Line 5 (heading) | Acronyms ONT and Illumina in the heading; ONT (Oxford Nanopore Technologies) is not expanded anywhere on this page. | Expand ONT as "Oxford Nanopore Technologies (ONT)" on first use in the body text. |
| low | ease-of-use | Line 6 | "long-read assemblies" and the short-read/long-read distinction assume familiarity; the value of polishing repeat regions is technical. | Add a brief plain-language note on why short-read polishing of a long-read assembly improves accuracy. |

### `docs/common_text/poppunk_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | accessibility | Line 6 | GPSC is expanded as "Global Pneumococcal Sequence Clusters" but GPS (used in 'GPS clusters', line 12, and links) is not expanded; k-mers (line 8) is jargon for a non-technical reader. | Expand GPS (Global Pneumococcal Sequencing) on first use and add a brief gloss for k-mers (short DNA subsequences of length k). |
| low | ease-of-use | Line 8 | "variable-length k-mers to distinguish between sample divergences in shared genomic content" and "pairwise distance distributions" are technical with no plain-language summary. | Add a one-line plain summary, e.g. "In short, PopPUNK groups genetically similar samples into clusters." |

### `docs/common_text/porechop_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | accessibility | Line 11 | ONT acronym used without expansion (Oxford Nanopore Technologies). | Expand ONT on first use: "ONT (Oxford Nanopore Technologies) data". |
| low | ease-of-use | Line 11 | "finding and removing adapters" assumes the reader knows what sequencing adapters are and why removing them matters. | Add a short clause explaining adapters are extra synthetic sequences added during library prep that should be removed before analysis. |

### `docs/common_text/prokka_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 10 | The sentence references the 'GFF3 file' and uses 'reduction of the GFF3' without explaining what a GFF3 file is, which a non-technical reader may not know. | Add a brief gloss, e.g. 'the GFF3 file (a standard text format that lists each annotated feature and its location)'. |
| low | ease-of-use | Line 10 | Hedged, vague phrasing ('is likely the most versatile') and a dense clause make the guidance unclear about when to use other formats. | Rephrase more directly, e.g. 'The GFF3 file is the most complete output because it contains all generated annotation information. Other formats are also produced for cases where a simpler subset of this information is more convenient.' |

### `docs/common_text/qc_check_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 8 | The output codes QC_PASS, QC_NA, and QC_ALERT are introduced but 'QC_NA' is only loosely explained ('could not proceed') with no example of what would cause that state, which a non-technical reader may find confusing. | Add a short example of when QC_NA appears (e.g. 'if a required metric was not generated, so the comparison could not be made'). |
| medium | ease-of-use | Line 17 (theiacov section) | This is a single very long, dense paragraph covering segment prefixes, general vs. segment-specific thresholds, and fallback behavior. It is hard to follow for the target audience. | Break into a short intro plus a bulleted list (one bullet for segment prefixes, one for the 'segment_' general prefix, one for how specific vs. general thresholds interact, one for the blank-cell rule). |
| low | grammar | Line 28 | 'copy and pasting' is grammatically awkward as a gerund pair. | Change to 'copying and pasting'. |
| low | accessibility | Line 20 (theiaprok\|theiaeuk section) | Acronym 'GAMBIT' (referred to as 'the GAMBIT module') is used without expansion or a one-line description of what it does. | On first use, briefly describe GAMBIT, e.g. 'the GAMBIT module (which predicts the sample's species)'. |

### `docs/common_text/qualimap_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 6 | Redundant repetition of 'various metrics' in the same sentence ('various metrics including ... and various metrics analyzed across the reference'). | Reword the end of the sentence, e.g. '...GC content, and other metrics measured across the reference.' |
| low | ease-of-use | Line 5-6 | Acronyms BAM, GC content, and 'next-generation sequencing' are used without expansion; non-technical readers may not know what a BAM file is. | Briefly gloss BAM on first use (e.g. 'BAM files, which store how sequencing reads align to a reference') and consider expanding NGS. |

### `docs/common_text/quasitools_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | accessibility | Line 6 | Acronyms 'WHO', 'HIVDR', and 'HyDRA' appear without expansion on first use. | Expand on first use, e.g. 'World Health Organization (WHO)' and note HyDRA is the drug-resistance module; HIVDR = HIV drug resistance. |
| low | ease-of-use | Line 8 | 'variant calling' is used without explanation, which may be unfamiliar to a non-technical reader. | Add a brief gloss, e.g. 'variant calling (identifying positions where the sample differs from the reference)'. |

### `docs/common_text/quast_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | Metrics 'number of contigs', 'length of the largest contig', and especially 'N50' are listed without definition; N50 in particular is unintuitive for non-technical readers. | Briefly define N50, e.g. 'N50 (a measure of assembly contiguity: the length such that half of the assembly is contained in contigs of that length or longer)', and gloss 'contig' as a contiguous assembled sequence. |

### `docs/common_text/racon_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 8 | 'consensus algorithm', 'de novo DNA assemblies', and 'polishing' are used without explanation; a non-technical reader may not know what assembly polishing accomplishes. | Add a one-line description of what polishing does (correcting errors in a draft assembly using the reads) and gloss 'de novo'. |

### `docs/common_text/rasusa_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 30 (read_qc_trim section) | Acronym 'HRRT' is introduced in parentheses with no expansion, and other tool names (Kraken2, MIDAS, Trimmomatic, fastp) appear without context for what they do. | Expand HRRT to 'host-read removal' fully (it is given as the expansion-before-acronym, but HRRT itself is never spelled out) and consider a brief note that Kraken2/MIDAS classify reads and Trimmomatic/fastp trim reads. |
| low | grammar | Line 17 | Missing comma after the introductory conditional clause: 'If more than one of `--bases`, `--frac`, or `--num` is supplied the task will fail'. | Add a comma: '...is supplied, the task will fail...'. |
| low | ease-of-use | Line 6 | 'subsample sequencing reads to a specified coverage' uses 'coverage' as core jargon without defining it for non-technical readers. | Add a brief gloss of coverage/depth, e.g. 'coverage (the average number of times each position in the genome is sequenced)'. |
| low | ease-of-use | Line 43 (ont section) | '5 million basepairs (0.7 Mb higher than the average bacterial genome length)' mixes units ('million basepairs' vs 'Mb') within one sentence, which can confuse readers. | Use consistent units, e.g. '5 Mb (about 0.7 Mb higher than the average bacterial genome length)'. |

### `docs/common_text/raven_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 7 | Dense, jargon-heavy single sentence ('overlap-layout-consensus based assembler', 'pile-o-grams', 'graph simplification method based on graph drawings', 'polishes unambiguous graph paths') that a non-technical reader will not understand. | Lead with a plain-language summary of what Raven does (builds a draft genome from long reads) and either drop or briefly gloss the highly technical internal-method terms. |
| low | grammar | Line 17 | 'difficult to traceback the source' uses 'traceback' as a verb; it should be two words as a verb. | Change to 'difficult to trace back the source'. |
| low | grammar | Line 17 | Awkward parallelism: 'fail with cryptic "segmentation fault" ... errors or by failing to output an output file' mixes a prepositional phrase with a gerund and repeats 'output an output'. | Reword, e.g. 'Raven may fail with cryptic "segmentation fault" (segfault) errors, or it may finish without producing an output file.' |

### `docs/common_text/read_decontaminate_wf.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 8 (theiaviral section) | Very dense paragraph describing three different ways to supply a host genome plus the mapping process; the input-field mechanics (is_genome, is_accession) are hard to parse inline. | Break the three host-input options into a bulleted list and separate the 'what the task does' steps from the 'how to provide input' details. |
| low | grammar | Line 8 (theiaviral section) | 'directly inputted' and 'inputted' are used; 'inputted' is nonstandard/awkward English. | Use 'input' or 'provided' as the past participle, e.g. 'directly provided' / 'reads provided as input'. |
| low | ease-of-use | Line 21 (read_qc_trim section) | 'comma-delimited string of expected sequence headers' and the requirement that they 'exactly match sequence headers in the input' assumes the reader knows what a FASTA sequence header is. | Add a brief explanation of what a sequence header is and give a short example of the expected_contaminants format. |

### `docs/common_text/read_qc_trim_illumina_wf.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | accessibility | Line 6 | Acronyms PE and SE are used ('the PE and SE versions') without expansion on first use. | Expand on first use: 'paired-end (PE) and single-end (SE) versions'. |

### `docs/common_text/read_qc_trim_ont_wf.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | accessibility | Line 5-6 | Acronym 'ONT' is used throughout the title and body without expansion. | Expand ONT to 'Oxford Nanopore Technologies (ONT)' on first use. |

### `docs/common_text/read_screen_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| high | ease-of-use | Line 30 (theiaviral table, max_genome_size row) | The rationale for `max_genome_size` reads 'A sample will fail the read screening if the estimated genome size is smaller than `max_genome_size`'. For a maximum threshold this is logically wrong — it should fail when the size is LARGER than the maximum. | Change 'smaller than' to 'larger than' for the max_genome_size row. |
| high | ease-of-use | Line 42 (theiacov table, max_genome_length row) / Line 54 (theiaprok) / Line 66 (theiaeuk) / Line 79 (theiaeukont) | Inconsistent variable naming between sections: the intro bullets and theiaviral table use `min_genome_size`/`max_genome_size`, while the theiacov/theiaprok/theiaeuk tables use `min_genome_length`/`max_genome_length`. Same concept, two different names, which is confusing. | Standardize on one term (e.g. `min_genome_length`/`max_genome_length`) across all sections, or note explicitly that the two names refer to the same input. |
| medium | ease-of-use | Line 27 (theiaviral table, min_reads row) and Line 29 (min_genome_size row) | Numeric inconsistencies undermine reader trust: line 29 says the rationale is the Hepatitis delta genome '(1,700 bp)' but the listed default `min_genome_size` is 1500; the theiacov table (line 41) uses 1700 for the equivalent value. Line 27 min_reads rationale text differs in form from the theiacov min_reads rationale despite describing the same calculation. | Reconcile the theiaviral default (1500 vs 1700) with the stated rationale, and align the min_reads rationale wording/values with the theiacov table. |
| medium | ease-of-use | Line 6 | Acronyms/jargon 'PE', 'SE', 'ONT', 'fastq-scan', 'Mash sketching', and 'coverage' appear without explanation; 'Mash sketching to estimate the genome size' is opaque to non-technical readers. | Expand PE/SE/ONT on first use and add a brief plain-language note on what Mash sketching does (rapidly estimates genome size from the reads). |
| low | grammar | Line 26 (theiaviral table, estimated_genome_length row) | 'This is an approximate of the median RNA virus length' — 'an approximate' is used as a noun, which is incorrect. | Change to 'This is an approximation of the median RNA virus length' or 'This is approximately the median RNA virus length'. |
| low | grammar | Line 10 | Missing terminal period: 'A sample will fail the read screening if there are fewer than `min_basepairs` basepairs' ends without punctuation. | Add a period at the end of the sentence. |
| low | grammar | Line 9 | Inconsistent reference 'the reads1 or read2 files' — 'reads1' vs 'read2' mismatch (one plural, one not) and inconsistent with the rest of the doc's 'read1/read2'. | Standardize to 'read1 or read2 files'. |
| low | grammar | Line 21 (theiaviral section) | 'Default values vary between the ONT and PE workflow' — 'workflow' should be plural to agree with 'between ... and'. | Change to 'between the ONT and PE workflows'. |
| low | ease-of-use | Line 6 | Very long opening sentence ending in a colon that then leads into a numbered list; the clause 'Samples are run through all threshold checks, regardless of failures, and the workflow will terminate...' is buried at the end. | Split into two or three shorter sentences, separating the description of the checks from the pass/fail termination behavior. |

### `docs/common_text/referenceseeker_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | accessibility | Line 17 / 'For ReferenceSeeker to identify a genome' paragraph | ANI is used as an acronym (in 'referenceseeker_ani_threshold') and again referenced on line 19 ('ranks highest according to ANI') but is never expanded for a non-technical reader. | Expand on first use, e.g. 'Average Nucleotide Identity (ANI)', and briefly explain it measures how similar two genomes are. |
| low | grammar | Line 17 | Extra space before the closing parenthesis: 'default >= 0.95 )' and trailing whitespace at end of line. | Remove the space so it reads 'default >= 0.95)'. |
| low | ease-of-use | Line 9 'Databases that can be used...' | GSURI is used without definition; a non-technical reader will not know what a GSURI is or where to paste it. | Expand to 'Google Storage URI (GSURI)' and clarify it is the file path string shown in the bullet list below. |

### `docs/common_text/resfinder_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 5 heading and line 10 'AMR Genotyping' | AMR and XDR appear in the heading and body without expansion on first use. A non-technical public-health reader may not know these acronyms. | Expand on first use: 'antimicrobial resistance (AMR)' and 'extensively drug-resistant (XDR)'. XDR is later expanded on line 31, which should be moved earlier or duplicated at first use in the heading/intro. |
| low | grammar | Lines 52-53 | Inconsistent count: criteria list (lines 38-47) defines 7 criteria including one for taxon and 6 drugs, but the output descriptions say 'all 6 drugs in the XDR definition above'. The drug count (6) is consistent, but earlier the criteria say 'all 7 criteria'. This may read as contradictory to a careful reader. | Clarify wording so the relationship between '7 criteria' (1 taxon + 6 drug resistances) and 'all 6 drugs' is explicit, e.g. 'resistance to all 6 drugs (criteria 2-7)'. |
| low | ease-of-use | Line 10 | 'BioNumerics' is referenced as a known reference point for thresholds but is not explained; readers may not know it is a commercial bioinformatics platform. | Add a brief gloss, e.g. 'BioNumerics (a commonly used commercial bioinformatics analysis platform)'. |
| low | ease-of-use | Line 15 | Cross-reference '[see above](#taxonomic-assignment)' assumes the reader is viewing this fragment inside a larger assembled page; within the fragment alone the anchor has no context, and GAMBIT is referenced without expansion. | Confirm the anchor resolves in the rendered page; consider briefly noting GAMBIT is the taxonomic identification tool used upstream. |

### `docs/common_text/seqsero2_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 6 | Awkward phrasing 'predict serotypes including contamination between serotypes' — 'including contamination' is ambiguous (contamination is not a serotype). | Rephrase, e.g. 'predict serotypes and can also flag potential contamination between serotypes.' |
| low | grammar | Line 8 | Run-on / comma splice: 'and any assembled alleles are used to predict the sample's serotype, and can predict potential contamination.' Subject of 'can predict' is unclear (alleles do not predict). | Rephrase, e.g. '...used to predict the sample's serotype and to flag potential contamination.' |
| low | ease-of-use | Line 8 | 'k-mer' and 'allele micro-assembly' are technical terms used without explanation for a non-technical audience. | Add a one-line gloss of k-mer (short DNA subsequence of length k) on first use. |

### `docs/common_text/seqsero2s_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 8 | Same ambiguous phrasing as seqsero2: 'predict serotypes including contamination between serotypes.' | Rephrase to 'predict serotypes and flag potential contamination between serotypes.' |
| low | grammar | Line 10 | Comma splice / unclear subject: 'any assembled alleles are used to predict the sample's serotype, and can predict potential contamination.' | Rephrase to '...used to predict the sample's serotype and to flag potential contamination.' |
| low | ease-of-use | Line 13 | LPS is used without expansion ('outermost oligosaccharides of LPS (O antigen)'). | Expand to 'lipopolysaccharide (LPS)'. |
| low | ease-of-use | Line 12-13 | PulseNet 2.0 / PN2.0 and the White-Kauffmann-Le Minor (WKL) Scheme are referenced as background but not explained; non-technical readers may not know what PulseNet is. | Add a brief gloss for PulseNet (the national molecular subtyping network for foodborne disease surveillance). |

### `docs/common_text/seroba_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 6 | 'k-mer based method' used without defining k-mer for a non-technical audience. | Add a brief gloss of k-mer on first use. |
| low | ease-of-use | Line 5 heading and line 6 | Multiple tool/database names (PneumoCaT, KMC, ARIBA) are introduced in dense technical prose without explaining what each does at a level useful to a non-bioinformatician. | Consider briefly noting the role of each (e.g. ARIBA = an assembly/variant-calling tool) or simplifying the method paragraph. |

### `docs/common_text/shared_variants_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | accessibility | Lines 20-24 (example table) | A complex example table with many abbreviated column headers (CHROM, POS, FTYPE, NT_POS, AA_POS, LOCUS_TAG, etc.) is shown without textual explanation of what the columns mean. | Add a short legend or sentence explaining the key columns so screen-reader users and non-technical readers can interpret the table. |
| medium | ease-of-use | Line 16 | Dense single paragraph explaining the 1/0 encoding; the meaning of '1' and '0' is important but buried in a long sentence. | Break into a short bulleted list: 1 = variant detected; 0 = variant not detected OR insufficient coverage. |

### `docs/common_text/shigapass_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 5 heading and line 6 | EIEC is expanded on line 6 ('enteroinvasive E. coli (EIEC)') but appears unexpanded in the heading on line 5; readers scanning the heading see the acronym first. | This is acceptable since it's expanded in the body, but consider expanding in the heading too for clarity. |

### `docs/common_text/shigatyper_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 6 | Double space between '_ipaB_' and 'which' ('reports on the presence of _ipaB_  which is suggestive'). | Remove the extra space and add a comma: '_ipaB_, which is suggestive'. |
| low | ease-of-use | Line 6 heading and body | EIEC appears in the heading unexpanded; it is expanded on line 6 in the body, which is fine, but minimap2 (line 8) is named without context. | Optionally note minimap2 is a sequence alignment/mapping tool for readers unfamiliar with it. |

### `docs/common_text/shigeifinder_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 9 | The acronym 'EIEC' is used in the heading and 'BLAST' and 'BWA' appear without explanation. While EIEC is expanded as 'enteroinvasive E. coli' in the body, BLAST and BWA are sequence-alignment tools that a non-technical reader will not recognize. | Briefly note that BLAST and BWA are sequence-alignment tools used to compare the sample against reference gene databases, or link to a glossary. |
| low | ease-of-use | Line 7 | The input variable described is 'call_shigeifinder_reads_input', but the surrounding text refers to the 'shigeifinder_reads' task. A reader may be unsure which exact toggle to set. The phrasing 'set the call_shigeifinder_reads_input to be true' reads awkwardly. | Rephrase to 'set call_shigeifinder_reads_input to true' and confirm the variable name matches the task name to avoid confusion. |

### `docs/common_text/sistr_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 (task body) | First paragraph is one dense block introducing several undefined concepts (cgMLST, centroid alleles, QC module) for a non-technical audience. | Expand 'cgMLST' (core genome multilocus sequence typing) on first use and consider splitting the paragraph; briefly explain what 'centroid/representative alleles' means in plain terms. |
| low | accessibility | Line 6 | Inline link text '[the QC section of the SISTR documentation]' is acceptable, but acronym 'QC' is not expanded on first use. | Expand 'QC' to 'quality control (QC)' on first use. |
| low | ease-of-use | Line 8 | 'MASH sketches' and 'MinHash algorithm' are advanced bioinformatics jargon used without explanation. | Add a one-line plain-language note that Mash/MinHash is a fast method for estimating genome similarity. |

### `docs/common_text/skani_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | grammar | Line 7 | 'It is magnitudes faster than BLAST-based methods' - 'magnitudes faster' is informal/imprecise; the standard phrasing is 'orders of magnitude faster'. | Change to 'It is orders of magnitude faster than BLAST-based methods'. |
| medium | ease-of-use | Line 7 | 'average nucleotide identity (ANI)' is expanded but its meaning/purpose is not explained for a non-technical reader; 'approximate mapping method without base-level alignment' is jargon. | Add a brief plain-language note that ANI measures how similar two genomes are (as a percentage), used here to find the closest reference. |
| low | grammar | Line 14 | Bullet item 2 lacks terminal punctuation while bullets 1 and 3 are full sentences with periods (inconsistent), and 'within the FASTA' is missing a noun (within the FASTA file). | Add a period and change to 'within the FASTA file' for consistency. |

### `docs/common_text/skesa_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 8 | 'de novo assembler', 'k-mer', 'haploid genomes', and 'repeat regions' are technical terms presented without explanation for a non-technical audience. | Add a brief plain-language note explaining that de novo assembly reconstructs a genome from sequencing reads without a reference. |

### `docs/common_text/snippy_core_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 8 | Acronyms 'SNPs' and 'INDELs' are used without expansion. For a non-technical audience these should be defined on first use. | Expand to 'SNPs (single-nucleotide polymorphisms) and INDELs (insertions and deletions)' on first use. |
| low | grammar | Line 14 | Double space between 'genomes' and 'genomes,' ... actually 'tuberculosis_  genomes' contains a double space before 'genomes'. | Remove the extra space: 'Mycobacterium tuberculosis genomes'. |
| low | ease-of-use | Line 14 | 'a bed file' is referenced without explanation of what a BED file is or how to make one, which a non-technical user masking regions would need. | Add a brief note describing what a BED file is (a plain-text file listing genome regions/coordinates) and where to learn how to create one. |
| low | ease-of-use | Lines 8, 11 | Phrase 'all the samples we'd like in our tree' uses a casual contraction and first-person plural inconsistent with the more formal tone of the rest of the documentation. | Rephrase to 'all the samples you want to include in the tree' for clarity and consistent voice. |

### `docs/common_text/snippy_qc_concatenation_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | accessibility | Line 8 | Acronyms 'QC' and 'TSV' (the latter appears in the related snippy_variants file) are not expanded; 'QC' should be defined on first use in this fragment. | Expand 'QC' to 'quality control (QC)' on first use. |
| low | ease-of-use | Line 29 | Single very long sentence (the tip block) packs multiple ideas (insights, monitoring, recommendation) together, making it dense for the target audience. | Split into two sentences: one stating what the metrics show, and one giving the recommendation to review them before phylogenetic analysis. |

### `docs/common_text/snippy_variants_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Lines 9, 12 | Acronyms 'SNPs', 'MNPs', and 'INDELs' are used without expansion. A non-technical reader will not know what MNPs are. | Expand on first use: 'SNPs (single-nucleotide polymorphisms), MNPs (multi-nucleotide polymorphisms), and INDELs (insertions and deletions)'. |
| low | grammar | Line 132 | Sentence states 'B9J08_000401 (which is FLO8)' but the search-term table at line 80 maps B9J08_000401 to FLO8, while the example output text on line 129 labels the same locus as 'hypothetical protein' - potentially confusing, and the example earlier says ERG11 missense at position 143 while the gene-name mapping is consistent. Verify the example interpretation matches the locus mapping to avoid contradiction. | Reconcile the example so the locus-to-gene-name mapping (FLO8 vs hypothetical protein) is consistent and clearly explained. |
| low | accessibility | Line 52 | 'AMRSearch' is referenced as a comparison without expansion or context; a reader unfamiliar with it cannot understand the contrast being drawn. | Briefly state what AMRSearch is (an antimicrobial resistance detection tool) when contrasting its behavior. |
| low | ease-of-use | Line 15 | Instruction about 'query_gene' requires the user to know what a GenBank file is and how the 'snippy_results' column is structured; assumes knowledge the audience may lack. | Briefly explain what a GenBank file is (an annotated reference sequence format) and give a short example of a valid query string. |

### `docs/common_text/snp_dists_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 9 | Very brief description uses 'multifasta alignment' and 'pairwise SNP distances' without explaining what these are or why the matrix is useful (no 'why/when to use this' guidance). | Add a sentence explaining that the SNP-distance matrix shows how many single-nucleotide differences separate each pair of samples, useful for assessing relatedness in outbreak investigations. |
| low | grammar | Line 9 | Trailing whitespace after the sentence (line ends with 'distances. ' plus a space). | Remove the trailing space. |

### `docs/common_text/snp_sites_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 11 | 'Gubbins' is referenced without any explanation of what it is or what 'recombination' means in this context, which assumes knowledge the audience may lack. | Add a brief note that Gubbins is a tool that detects recombinant regions, and that excluding them improves phylogenetic accuracy. |
| low | ease-of-use | Line 5 (heading) and Line 14 | The task heading is 'Genome Filtering' but the description says SNP-sites 'identifies variants' and returns 'only the sites with SNPs' - the relationship between 'filtering' and 'extracting SNP sites' may be unclear to non-technical readers. | Add a one-line purpose statement clarifying that SNP-sites filters an alignment down to only the variable (SNP) positions to reduce data size for phylogenetics. |

### `docs/common_text/sonneityper_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | Acronym 'ONT' appears in the heading without expansion, and gene/mutation notations (gyrA S83L, parC S80I) are highly technical with no plain-language context about what they mean for the user. | Expand 'ONT' to 'Oxford Nanopore Technologies (ONT)' on first use, and briefly note that the listed mutations are known markers of quinolone (antibiotic) resistance. |
| low | grammar | Line 8 | British spelling 'analyses' used as a verb ('a tool, Mykrobe, that analyses'). The repo standard is American English, where the verb is spelled 'analyzes'. | Change 'analyses' to 'analyzes'. |

### `docs/common_text/spades_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 20 | 'MetaviralSPAdes pipeline consists of three independent steps' is missing the article 'The' at the start of the clause. | Change to 'The MetaviralSPAdes pipeline consists of three independent steps'. |
| low | grammar | Line 20 | Run-on sentence listing the three steps with inconsistent punctuation (commas separating steps, 'and' before the last) makes the dense sentence hard to parse for the target audience. | Convert the three steps into a bulleted list (ViralAssembly, ViralVerify, ViralComplete) for readability. |
| low | accessibility | Line 34 | The techdetails table is titled 'MetaviralSPAdes Technical Details', but in TheiaProk the task runs plain SPAdes (not MetaviralSPAdes). The fixed title may mislead TheiaProk readers. | Use a more general title such as 'SPAdes Technical Details' so it is accurate for both TheiaProk and TheiaViral contexts. |
| low | ease-of-use | Line 13 | 'de Bruijn graphs', 'de novo assembly', and 'Illumina short reads' are technical terms with no plain-language explanation for a non-technical audience. | Add a brief note that de novo assembly reconstructs a genome from sequencing reads without using a reference. |
| low | ease-of-use | Line 22 | 'CheckV quality assessment' is mentioned without explaining what CheckV is. | Briefly note that CheckV is a tool that estimates the completeness and quality of viral genome assemblies. |

### `docs/common_text/spatyper_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | accessibility | Line 6 | The text warns the Ridom SpaServer link 'is typically considered unsafe by modern web browsers' but provides no guidance on what the user should do, which could cause confusion or hesitation. | Add a brief reassurance/instruction (e.g., that the warning is due to the site's certificate and the tool uses the database internally, so users typically do not need to visit the link). |
| low | ease-of-use | Line 8 | Acronym 'MRSA' is used and expanded inline, which is good; however 'spa-type' and 'protein A gene' concepts assume some background. This is minor but a one-line 'why' (what a spa-type tells you) would help. | Optionally add a short note that a spa-type is a strain identifier used to compare isolates during outbreak tracking. |

### `docs/common_text/srst2_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 5 (task title) / Line 6 | SRST2 is never expanded. The acronym appears in the title and body with no definition for a non-technical reader. | Expand on first use, e.g. 'SRST2 (Short Read Sequence Typing for Bacterial Pathogens)'. |
| low | grammar | Line 21 | Subject-verb agreement: 'The presence or absence of specific genes are used' should be 'is used' (subject is 'presence or absence'). | Change 'are used' to 'is used'. |
| low | accessibility | Line 8 toggle 'Resistance Gene Database' / table lines 9-17 | The toggle is labeled 'Resistance Gene Database' but the table contains species markers, toxin, biotype and serogroup genes, not resistance genes. Mislabeling can confuse and mislead, especially for screen-reader users navigating by the toggle label. | Rename the toggle to something accurate, e.g. 'Vibrio Characterization Gene Database'. |
| low | ease-of-use | Line 6 | 'target sequences that are traditionally used in PCR methods' and the concept of a 'database of target sequences' may be unclear to non-technical readers without a brief 'why'. | Add a short clause explaining the purpose, e.g. 'a curated set of genes used to characterize Vibrio (such as confirming species and toxin production)'. |
| low | ease-of-use | Line 21 | Jargon like 'biotype' and 'serogroup' is used without definition for a non-technical audience. | Add a one-line gloss for biotype and serogroup, or link to a glossary. |

### `docs/common_text/staphopia_sccmec_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | MRSA is used without expansion. A non-technical reader may not know the acronym, and the sentence implies mecA causes MRSA without stating what MRSA stands for. | Expand on first use: 'methicillin-resistant Staphylococcus aureus (MRSA)'. |
| low | grammar | Line 6 | Redundant 'also'/'as well': 'has also been found to confer resistance to non-beta-lactam drugs as well'. | Remove one redundancy, e.g. 'has also been found to confer resistance to non-beta-lactam drugs'. |
| low | grammar | Line 6 | Slight imprecision: 'a methicillin-resistant mecA gene'. The gene confers methicillin resistance; the gene itself is not 'methicillin-resistant'. | Reword to 'a methicillin-resistance gene (mecA)' or 'the mecA gene, which confers methicillin resistance'. |
| low | ease-of-use | Line 8 | Hamming Distance is linked but not briefly defined; a non-technical reader gets no inline sense of what it measures here. | Add a short gloss, e.g. 'the Hamming Distance (a count of how many positions differ between two sequences)'. |

### `docs/common_text/stxtyper_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | accessibility | Line 8 | Vague/awkward link text: 'Please see this review paper that describes shiga toxins in great detail' — the whole sentence is the link, starting with 'Please see this'. Better link text should describe the destination concisely. | Make the descriptive title the link text, e.g. 'See [this review of Shiga toxins](...) for more detail.' |
| medium | ease-of-use | Line 6 / Line 8 | STEC is introduced on line 8 (E. coli appears earlier on line 6). E. coli, Shigella, Klebsiella are italicized inconsistently; also 'E. coli' on line 6 is not italicized while species names elsewhere are. | Italicize species names consistently (E. coli, Shigella, Klebsiella) and ensure STEC is expanded at first use (it is expanded as Shiga-toxin-producing E. coli). |
| low | grammar | Line 6 | Inconsistent capitalization of 'shiga toxin' vs 'Shiga toxin'. Title and line 5 use 'Shiga'; line 6 uses lowercase 'shiga toxin' and 'shiga-toxin-producing'. | Capitalize 'Shiga' consistently throughout (Shiga toxin, Shiga-toxin-producing). |
| low | grammar | Line 13 | Inconsistent tool name casing: 'Stxtyper' on line 13 vs 'StxTyper' elsewhere; also a double space after 'composition.' | Standardize to 'StxTyper' and remove the double space. |
| low | accessibility | Line 15 | Link text 'this output format.' is vague and includes the period inside the link, which is poor for screen readers. | Use descriptive link text such as '[the StxTyper output format documentation]' and place the period outside the link. |
| low | ease-of-use | Line 13 | Jargon 'operon', 'subunits', 'amino acid composition' used without definition for non-technical readers. | Add a brief gloss for 'operon' (a group of genes transcribed together) or link to a definition. |

### `docs/common_text/taxon_table_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | grammar | Line 9 and Line 30 | 'This task is incompatible when running TheiaViral_Panel/TheiaProk on the command-line' is awkward phrasing — a task is not 'incompatible when running'. | Reword, e.g. 'This task is not compatible with command-line runs of TheiaViral_Panel; it is designed specifically for Terra.' |
| low | grammar | Lines 37, 41, 61 | Inconsistent plural of 'taxon': 'samples of each taxa', 'samples matching a taxa'. 'Taxa' is plural; with 'each' / 'a' the singular 'taxon' is correct. | Use 'each taxon' and 'matching a taxon'. |
| low | accessibility | Lines 8, 15, 29, 36, 60 | Multiple admonitions use empty titles (`!!! tip ""`). Empty titles provide no context for screen-reader users navigating by landmark/heading. | Provide concise titles for each tip admonition describing its content. |
| low | ease-of-use | Line 18 / Line 41 | 'leftmost column' and 'rightmost column' describe a two-column table; minor, but pairing the description with the header names ('taxon' and 'taxon_table') would be clearer. | Reference the actual column headers, e.g. 'the taxon name in the `taxon` column and the destination data table name in the `taxon_table` column.' |

### `docs/common_text/tbp_parser_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 5 / Line 8 | LIMS is used (line 10) without expansion; also TBProfiler/JSON appear in the title aimed at a non-technical audience without context for what 'parsing JSON outputs' means. | Expand LIMS on first use ('Laboratory Information Management System (LIMS)') and briefly state in plain terms what tbp-parser does for the user. |
| low | grammar | Lines 8, 13 | Trailing whitespace at end of lines 8 and 13 (cosmetic). | Remove trailing spaces. |
| low | accessibility | Line 8 | Link text 'please examine the documentation for this tool' is a vague directive rather than describing the destination. | Use descriptive link text such as '[the tbp-parser documentation]'. |
| low | ease-of-use | Line 12 | ONT acronym used without expansion ('Oxford Nanopore Technologies') for a non-technical reader. | Expand ONT on first use or link to a definition. |

### `docs/common_text/tbprofiler_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 10 | Trailing whitespace at end of line 10 (cosmetic). | Remove trailing spaces. |
| low | ease-of-use | Line 8 | Heavy unexplained jargon for a non-technical audience: BWA, mem, Minimap2, GATK, BCFtools, FreeBayes, LoFreq, Pilon, 'calls variants'. No 'why' or plain-language summary. | Add a brief plain-language sentence describing the overall process (aligns reads to a reference, finds mutations, compares to a known database) before the tool names. |
| low | ease-of-use | Line 10 | Dense paragraph mixing format guidance and capability list; 'not very human readable' is informal/awkward. | Split into two sentences/bullets and reword to 'the JSON file is difficult for people to read directly'. |
| low | ease-of-use | Line 6 | 'sub-lineages' and 'drug resistance-associated mutations' assume domain knowledge; acceptable but a one-line 'what TBProfiler is for' would help the target audience. | Add a short 'when to use this' note (e.g., for TB drug-resistance surveillance). |

### `docs/common_text/theiaviral_inputs.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | grammar | Line 9 | Redundant/awkward wording: 'the first retrieved genome corresponding to that taxon is retrieved' — 'retrieved' used twice in one clause. | Reword, e.g. 'the first genome corresponding to that taxon is used.' |
| medium | ease-of-use | Line 9 | Very dense single sentence packing taxon/FASTA/accession options plus four coupled boolean variables and PE/ONT distinctions; hard to follow for non-technical readers. | Break into a short list separating the three input types and which boolean each requires. |
| medium | ease-of-use | Line 14 | HRRT acronym used without expansion ('Human Read Removal Tool'); also Skani, RefSeq, de novo assembly used without brief context for the target audience. | Expand HRRT on first use and add brief glosses or links for Skani and RefSeq. |
| low | grammar | Line 22 | Missing terminal period at the end of the last sentence ('...Kraken2 for Illumina and Metabuli for ONT)'). | Add a closing period. |
| low | grammar | Lines 14, 15 | Trailing whitespace at end of lines 14 and 15 (cosmetic). | Remove trailing spaces. |
| low | ease-of-use | Line 22 | Single very long, dense paragraph covering read_extraction_rank with multiple nested caveats; difficult for the audience to parse. | Break into shorter sentences or bullets (default behavior, why raise the rank, the risk of raising it too far, how to refine). |
| low | ease-of-use | Line 12 / Line 19 | 'allele', 'allelic variant', 'consensus assembly', 'mapping quality', 'depth' used without definition; non-technical readers tuning min_allele_freq/min_depth/min_map_quality may not understand the impact. | Add brief plain-language definitions or a short example of what changing each parameter does to results. |

### `docs/common_text/trimmomatic_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | grammar | Line 8 | Parameter name mismatch likely a typo: the window size is 'specified with `trim_window_size`' but the quality threshold is 'trimmomatic_window_quality'. The inconsistent prefix (trim_ vs trimmomatic_) reads like an error. | Verify and standardize the parameter name (likely `trimmomatic_window_size`) so prefixes are consistent. |
| low | grammar | Line 20 | Missing word: 'using the `trimmomatic_adapter_fasta` and `trimmomatic_adapter_trim_args` respectively' — should be '...parameters respectively'. | Add 'parameters' after the second variable name. |
| low | ease-of-use | Line 8 | 'sliding window', 'average quality', 'seed mismatches', 'palindrome clip threshold', 'simple clip threshold' (lines 15-17) are advanced terms with no plain-language explanation for the audience. | Add a one-line plain description of what a sliding window quality trim does, and note these ILLUMINACLIP values are advanced settings most users can leave at default. |
| low | ease-of-use | Line 20 | 'colon-delimited values that come after the adapter fasta file in the ILLUMINACLIP argument' assumes the reader can parse the ILLUMINACLIP syntax. | Reference the labeled example in the code block above (seed mismatches : palindrome clip : simple clip) to anchor the explanation. |

### `docs/common_text/ts_mlst_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 35 | Double space after 'genera,' ('the _Escherichia_ and _Leptospira_ genera,  _Acinetobacter'). | Replace the double space with a single space. |
| low | grammar | Line 45 | The Original Publication(s) table cell is missing a closing pipe/period — the row ends with the link and no terminal '\|' on the visible line, and the sentence lacks a period. | Ensure the table row is properly closed with a trailing pipe and add a closing period if intended. |
| low | ease-of-use | Line 10 | Dense, jargon-heavy paragraph (loci, allele, allelic profile, sequence type, clonal complexes, recombinational events) for a non-technical audience. | Consider a short bulleted breakdown of the key terms (locus, allele, ST, clonal complex) for readability. |
| low | ease-of-use | Line 13 | 'some MLST schemes are insufficiently robust' is vague — it states a limitation without telling the user what to do about it. | Add brief guidance, e.g. recommend verifying ambiguous results or consulting PubMLST for the scheme's status. |

### `docs/common_text/vadr_flu_segments_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 5 / Line 7 | VADR is not expanded in this fragment and 'multifasta', 'contigs', '.tar.gz', '.sqc file' are technical terms with no plain-language context for the audience. | Expand VADR on first use and add brief glosses (e.g., 'multifasta — a single file holding multiple sequences'). |
| low | accessibility | Line 9 ('Note:') | Important caveat about unreliable results is in a plain 'Note:' paragraph rather than an admonition, making it easy to overlook and less distinguishable for screen-reader navigation, unlike the admonition style used elsewhere. | Convert to a '!!! warning' or '!!! note' admonition for consistency and emphasis. |
| low | ease-of-use | Line 5 / Line 7 | No 'why/when to use this' — the fragment explains what it does but not why a user would run it (separating flu segments). | Add a one-line purpose statement (e.g., 'Use this when you need each influenza segment as a separate file, for example for per-segment analysis or submission'). |

### `docs/common_text/vadr_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | accessibility | Line 6 | Link text 'organism-specific parameters and logic section' is descriptive (good), but the surrounding instruction relies on the reader knowing what 'default models/parameters' means with no context. | Briefly note what the linked section covers (which models VADR uses per organism). |
| low | ease-of-use | Line 8 | 'Fatal alerts' vs 'non-fatal alerts' and 'designated as passing sequences' could be clearer; the sentence 'non-fatal alerts are designated as passing sequences' conflates alerts with sequences. | Reword to 'sequences with only non-fatal alerts still pass, but may require further investigation.' |
| low | ease-of-use | Line 6 | GenBank is referenced without a one-line description for non-technical readers. | Add a brief gloss, e.g. 'GenBank (NCBI's public sequence database)'. |

### `docs/common_text/versioning_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 6 | The phrase 'captures the workflow version from the GitHub (code repository) version' is awkward and repeats 'version' twice, making it hard to parse for a non-technical reader. | Rephrase to something like: 'The versioning task records which version of the workflow was used, based on the version tag in the GitHub code repository.' |

### `docs/common_text/vibecheck_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 6 | The acronyms 'SRST2' and 'SNPs' and the term 'variant frequency demixing' are used without definition. A non-technical public health scientist would not know these terms. | Expand 'SNPs' as 'single nucleotide polymorphisms (SNPs)' on first use, briefly note what SRST2 does (e.g. 'the SRST2 typing task'), and add a short plain-language gloss for 'variant frequency demixing'. |
| low | grammar | Line 8 | The sentence has a doubled construction: 'estimating lineage abundances using Freyja by using a database' — 'using ... by using' is awkward. | Revise to: 'estimating lineage abundances using Freyja with a database built from canonical SNPs that define the known lineages.' |
| low | ease-of-use | Line 5 (task title) | The 'O1' designation in the title is jargon that is not explained until later, and even then only briefly. A reader unfamiliar with V. cholerae serogroups will not know what O1 means. | Add a brief parenthetical defining O1 as the serogroup designation (e.g. 'O1 serogroup') on first use. |

### `docs/common_text/virulencefinder_task.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | grammar | Line 6 | Broken/run-on sentence. 'VirulenceFinder, from the Centre for Genomic Epidemiology (CGE) in TheiaProk is only run on assembly files...' The clause 'in TheiaProk' is misplaced and there is no comma closing the parenthetical phrase, making the sentence ungrammatical. | Restructure, e.g.: 'VirulenceFinder, from the Centre for Genomic Epidemiology (CGE), is run within TheiaProk only on assembly files, due to discordant results observed when using read files on the web application versus the command line.' |
| low | ease-of-use | Line 6 | Acronyms 'BLAST' and 'KMA' are used without expansion or explanation for a non-technical audience. | Briefly note that these are sequence-comparison tools (e.g. 'uses the alignment tools BLAST and KMA') so readers understand their role. |

### `docs/contributing/code_contribution.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | accessibility | Code block lines 47-58 (Bracket and Spacing Conventions) | The 'Correct' example block has inconsistent indentation (the opening `input {` is indented two spaces less than the comment above it), which could confuse readers relying on the example for formatting guidance. | Align the indentation of the 'Correct' example so the comment and code start at the same column. |
| low | ease-of-use | Line 41 | Nested parentheses make this instruction hard to parse: 'Use a single space when defining variables (`this = that` _not_ `this= that` (unless a bash variable where `this=that` is required))'. | Split into two sentences or restructure to avoid nested parentheses, e.g.: 'Use a single space around the equals sign when defining variables (`this = that`, not `this= that`). The exception is bash variables, which require no spaces (`this=that`).' |
| low | ease-of-use | Line 123 | The instruction 'Expose Docker as an input, an output (if versioning information not available), and runtime variable' is terse and the parenthetical condition is ambiguous about which item it applies to. | Clarify which element the 'if versioning information not available' condition modifies, e.g. 'expose Docker as an input and a runtime variable, and also as an output if version information is not otherwise captured.' |

### `docs/contributing/doc_contribution.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Lines 184-185 ('Adding a Page for a New Workflow') | Steps reference editing the `mkdocs.yml` file, but the intro (line 3) states the project now uses Zensical rather than MkDocs. This is inconsistent terminology that could confuse a contributor about which configuration file to edit. | Verify and update the filename if the navigation config has changed under Zensical, or add a note clarifying that `mkdocs.yml` is still the active config file despite the move to Zensical. |
| low | grammar | Line 27 | Missing word: 'allow you preview changes' should be 'allow you to preview changes'. | Change to 'allow you to preview changes'. |
| low | grammar | Line 174 | Missing article: 'Add a page with the title of the workflow to appropriate subdirectory' is missing 'the' before 'appropriate'. | Change to 'to the appropriate subdirectory'. |

### `docs/getting_started/commandline.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 4 ('What is WDL?') | WDL is expanded inconsistently: here it is 'Workflow Development Language', but in Step 3 (line 106) it is 'Workflow Description Language'. Inconsistent terminology for the same acronym confuses readers. | Use one consistent expansion throughout. The official expansion is 'Workflow Description Language'. |
| low | accessibility | Line 174 and line 220 (headings inside toggles) | Heading-level jump: the page uses H2 (##) for steps, but inside the toggle admonitions H5 (#####) headings appear (e.g. 'Tips for monitoring workflow progress', 'Canceling a Running Workflow'), skipping H3 and H4 levels. | Use a heading level consistent with the surrounding hierarchy (e.g. H3 or H4) to avoid skipped levels for screen-reader navigation. |
| low | ease-of-use | Line 4 | 'Frank has put together a great video' references a person ('Frank') by first name with no context; a new external reader will not know who Frank is. | Provide context, e.g. 'A member of our team has put together a great video' or identify Frank's role. |
| low | ease-of-use | Heading 'Step 2: Install docker and miniWDL' (line 73) | Inconsistent capitalization/spelling of tool names: 'docker' vs 'Docker' and 'miniWDL' vs 'miniwdl' vs 'miniWDL' are used inconsistently across the page (e.g. line 73 'docker and miniWDL', line 75 'Docker and miniwdl'). | Standardize to 'Docker' and 'miniwdl' (the tool's canonical spellings) throughout for consistency. |

### `docs/getting_started/terra.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 51 ('TheiaX workflows') | The placeholder term 'TheiaX' is used as an umbrella name without first explaining that it stands for the family of Theia* workflows (TheiaCoV, TheiaProk, etc.). A new reader may not realize 'TheiaX' is a collective shorthand. | Add a brief note on first use, e.g. 'The TheiaX workflows (a collective name for TheiaCoV, TheiaProk, TheiaEuk, and TheiaMeta) ...'. |
| low | ease-of-use | Line 4 ('Our Approach') | Acronyms used without expansion for a non-technical audience: 'INSDC' (line 38), 'SRA', 'ENA', 'DRA' (line 38), and 'QC' (line 63) appear without definition on first use. | Expand acronyms on first use, e.g. 'International Nucleotide Sequence Database Collaboration (INSDC)' and 'quality control (QC)'. The repository names already link out, but a brief gloss would help. |

### `docs/guides/custom_organisms.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 1 (title) / overall | The page title is 'Guide to Running Custom Organisms on TheiaCoV' but the opening lines and the 'Workflow Recommendations' section steer users toward TheiaViral instead, calling the TheiaCoV content 'legacy information.' The purpose and 'when to use this' guidance is muddled — a reader is not clearly told when they should follow this guide versus going to TheiaViral. | Add a short upfront framing box stating who this guide is for and when to use TheiaCoV custom-virus support versus TheiaViral, so non-technical users can decide quickly. |
| low | ease-of-use | Line 32 vs line 19-29 | Inconsistent count of supported organisms: the bulleted list (lines 22-29) shows eight organisms but line 32 says 'currently support eight organisms (see above)' — counting the list yields 8 entries, but this should be verified. The redundant restating ('These workflows currently support the following organisms' then 'These workflows currently support eight organisms') is also repetitive. | Remove the redundant sentence on line 32 or merge it with the list intro, and double-check the organism count matches the list. |
| low | ease-of-use | Line 44 (paragraph beginning 'All non-influenza default organisms') | Very dense single paragraph packed with many tool names (HRRT, Kraken2, trimmomatic, fastp, bbduk, fastq_scan, FastQC, bwa, iVar, samtools) and steps, which is hard for a non-technical reader to follow. | Break into a short numbered list of pipeline steps (read removal, taxonomic profiling, trimming, adapter removal, QC, mapping, primer trimming, variant calling, consensus) to improve readability. |
| low | ease-of-use | Line 34 and line 50 | References to 'Figure 2 below' (line 34) and figure numbering rely on the reader tracking figure numbers; the assembly flowchart is labeled figure-2 but the first figure is figure-1. For a non-technical reader, the cross-reference 'see Figure 2 below' without a hyperlink may be hard to locate. | Make 'Figure 2' a clickable cross-reference/anchor link, consistent with how figures are referenced elsewhere (e.g. the GAMBIT page links to #figure2). |

### `docs/guides/gambit.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | accessibility | Line 30 (card) and lines 477-547 (GAMBIT Fungal Databases sections) | The landing card links to a 'Latest Database Version' including 'GAMBIT Fungal Database v1.0.0', and full fungal database sections remain, but the repository's recent commit history indicates the fungal content was being removed. If fungal support was intentionally removed elsewhere, these references are now orphaned/inconsistent. Worth verifying the fungal sections are still intended to be present. | Verify whether the GAMBIT Fungal Database sections and links should remain; if fungal content was removed elsewhere, update or remove these references for consistency. |
| low | accessibility | Line 584 / line 586 (Figure 1 description) | The figure description relies on color to convey meaning ('The green histogram...', 'The blue histogram...', 'The red histogram...', 'dashed blue line') without a non-color textual cue for which distribution is which. Screen-reader and color-blind users cannot map colors to meaning. | Add a non-color identifier for each distribution (e.g. label them A/B/C or by position) in addition to color, or ensure the alt text conveys which distribution represents intra-species vs inter-species without relying solely on color. |
| low | accessibility | Lines 695, 740, 772 (DBeaver instructions) | Instructions reference 'click the orange arrow to execute SQL statements' — color-only reference to a UI control with no alternative descriptor. | Add a non-color descriptor, e.g. 'the orange (Execute) arrow button' or describe its location/icon, so the instruction does not rely on color alone. |
| low | ease-of-use | Line 598 ('GAMBIT distances correlate with sequence identity') | Acronym 'ANI' is expanded as 'Average Nucleotide Identity' here, which is good, but it is used earlier without expansion in some contexts; also 'Spearman correlation' (line 610) and the negative correlation values are presented without a plain-language explanation of what a value near -0.97 means for a non-technical reader. | Add a brief plain-language note that a Spearman value close to -1 indicates a very strong inverse relationship (higher ANI means lower GAMBIT distance). |
| low | ease-of-use | Line 576 ('Built-in Species Distance Threshold') | Sentence fragment / grammar: 'Some species are not well separated from their closest sister taxon and, in some cases, even overlap. Such as the case of Escherichia coli and Shigella sonnei...' — 'Such as the case of...' is a sentence fragment. | Join to the prior sentence, e.g. '...even overlap, such as Escherichia coli and Shigella sonnei in GAMBIT's Prokaryotic Database.' |
| low | ease-of-use | Lines 22-32 (cards) and line 18 | Acronyms 'GTDB' and 'GS URI' are used in the cards and database sections before being expanded (GTDB is first spelled out only at line 294). Non-technical readers encounter 'GTDB Database' early without knowing it means Genome Taxonomy Database. | Expand 'GTDB' as 'Genome Taxonomy Database (GTDB)' on its first appearance, and briefly note that a GS URI is a Google Cloud Storage path. |

### `docs/guides/gambit_database.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 41 | Double space inside the sentence: 'genomes sourced from RefSeq and GenBank. Automation is available through the GAMBITdb which takes...' — also 'through the GAMBITdb' reads awkwardly (GAMBITdb is a tool name, not 'the GAMBITdb'). | Change to 'Automation is available through GAMBITdb, which takes as input...' (remove 'the'). |
| low | grammar | Line 55 | Awkward passive phrasing: 'Species missing from previous iterations of the database were attempted to be added with looser quality criteria.' | Rephrase to active voice, e.g. 'We attempted to add species missing from previous database iterations using looser quality criteria.' |
| low | grammar | Line 63 | Comma splice / misplaced 'however': 'GTDB proposes new genus/species based on ANI however, these are not currently formally accepted by the community.' | Repunctuate: 'GTDB proposes new genera/species based on ANI; however, these are not currently formally accepted by the community.' |
| low | accessibility | Line 3 / Line 1 | Acronyms 'GTDB', 'ANI', 'N50', and 'GS URI' appear without expansion. GTDB is linked but never spelled out on this page; ANI (line 63) and N50 (line 135) are unexpanded jargon for a non-technical audience. | Expand on first use: 'Genome Taxonomy Database (GTDB)', 'Average Nucleotide Identity (ANI)', and add a brief note that N50 is a genome-assembly contiguity metric. |
| low | accessibility | Line 155 (GAMBIT Metadata Database Schema figure) | The figure caption text below the image reads 'Scheme for GAMBIT's metadata database' — 'Scheme' should be 'Schema' (typo) and the alt text, while descriptive, could note this is a reference diagram for understanding the .gdb file structure. | Correct 'Scheme' to 'Schema' in the figure caption on line 157. |
| low | ease-of-use | Line 47 (Curation of Genomic Data) | Very long, dense sentence: 'This is because many genomes are identical or nearly identical to others in the public archives and there is no benefit to having multiple examples, for example, in an outbreak dozens of identical genomes could be sequenced and after the first, there is no new information for classification.' This run-on with two 'for example' style clauses is hard to follow. | Split into two sentences and remove the redundancy, e.g.: 'Many genomes in public archives are identical or nearly identical, so storing multiple copies adds no value. For example, in an outbreak dozens of identical genomes may be sequenced, but after the first, none add new information for classification.' |

### `docs/guides/phylogenetics.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| high | ease-of-use | Line 45 ('heterologous sites'); also 'HGT' on Line 61 | Undefined jargon for a non-technical audience. 'heterologous sites' is highly technical and undefined; 'HGT' is an unexpanded acronym on first use. | Define or simplify 'heterologous sites' (e.g., 'sites where the sequence differs from the reference'). Expand HGT to 'horizontal gene transfer (HGT)' on first use. |
| medium | grammar | Line 46 | Missing space before the en/hyphen dash: 'should be closely related- the same lineage or sequence type'. The hyphen is jammed against the preceding word and used where a colon or em dash is intended. | Change to 'should be closely related: the same lineage or sequence type' (or use a spaced em dash). |
| medium | accessibility | Line 5 (heading) and Line 39, 74 etc. | The H2 heading 'Broadly, there are two phylogenetic analysis methods' is a full sentence used as a heading, which is unusual for navigation and screen-reader heading lists. Headings should be concise labels. | Reword to a noun-phrase heading such as 'Two Main Phylogenetic Analysis Methods' and move the 'Broadly,' framing into body text. |
| medium | accessibility | Lines 58-65 (comparison table) and lines 97-111 (visualization table) | Large multi-column comparison tables are presented with no surrounding textual summary. Screen-reader users and non-technical readers get no plain-language overview of what the tables convey or how to choose. | Add a short introductory sentence before each table summarizing the key takeaway (e.g., 'Use this table to compare speed, accuracy, and input requirements across the workflows; in general, choose X for Y'). |
| medium | ease-of-use | Line 61 and throughout the comparison table (lines 61-65) | Several undefined acronyms/terms for non-technical readers: 'NJ tree', 'ML tree', 'kmer'/'k-mer', 'CDS', 'GFF3', 'bed file', 'SNP'. These appear without expansion. | Add a short glossary or expand on first use: NJ (neighbor-joining), ML (maximum likelihood), k-mer, CDS (coding sequence), and define SNP (single nucleotide polymorphism) the first time it appears. |
| medium | ease-of-use | Line 76 (end of paragraph before bullets) | Broken/incomplete sentence: 'Typically, SNP distance thresholds can' ends with no terminal punctuation and flows into a bulleted list, but reads as a sentence fragment. | End the lead-in clearly, e.g., 'Typically, SNP distance thresholds can be used to:' followed by the list. |
| low | typo | Line 84 (SNP threshold bullet list, lines 82-85) | The bulleted list under 'It can be difficult to determine SNP thresholds because of:' is indented with leading spaces under a paragraph, which in Markdown will likely render as a code block or fail to render as a proper list (the colon line is not a list item and the bullets are over-indented). | Remove the extra indentation so the dashes start at the left margin (or are nested correctly) to render as a real bulleted list. |
| low | grammar | Line 61, table cell for Mashtree | Missing period and abbreviation style: '(recombination, HGT, etc)' lacks a period after 'etc'. | Use 'etc.' with a period. |
| low | grammar | Line 87 | Uses curly/typographic apostrophes ('don’t', 'weren’t') while the rest of the document uses straight apostrophes elsewhere; also 'allowing inference of start of outbreak' is missing an article. | Make apostrophe style consistent, and change to 'allowing inference of the start of an outbreak'. |
| low | grammar | Line 101 (visualization table, MicrobeTrace 'Ease of use' cell) | Doubled space: 'change the visualization from network  to tree'. | Remove the extra space. |
| low | ease-of-use | Lines 50-56 ('Workflow recommendations') | Workflow names like 'Augur_Prep & Augur', 'kSNP3', 'Snippy_Streamline' are listed without linking to their dedicated pages, forcing a non-technical reader to search for them. | Link each workflow name to its documentation page so readers can quickly learn what each does. |

### `docs/guides/sops.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 1-3 (intro) and all section headings | The acronyms SOP and PHB are used without expansion in headings/tables; SOP is expanded once in the title but column references like 'PHB version compatibility' assume the reader knows PHB. | Add a one-line note in the intro defining PHB (Public Health Bioinformatics) and reaffirming SOP (Standard Operating Procedure) for first-time readers. |
| low | ease-of-use | Lines 61-65 ('Miscellaneous SOPs') | The 'Miscellaneous SOPs' section filters on categories including 'Data Import', which is also its own dedicated section (lines 21-27), creating overlap that may confuse readers about where to find a given SOP. | Clarify in a short sentence what 'Miscellaneous' contains, or remove the overlapping 'Data Import' filter from the miscellaneous section. |

### `docs/guides/versions.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 1-3 (intro) | Assumes the reader knows what 'Docker images', 'databases', and 'reference files' mean in this context; non-technical public health scientists may not understand why these defaults matter or when they would change them. | Add one or two plain-language sentences explaining what a Docker image is (a packaged version of a software tool) and noting that most users can leave these defaults unchanged. |

### `docs/index.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 7 | Says 'Workflows are available for viruses, bacteria, and fungi' but the repository's terminology elsewhere uses 'Mycotics' for fungi (e.g., workflows_kingdom.md). Inconsistent terminology for the same concept could confuse readers cross-referencing pages. | Use consistent terminology, e.g., 'viruses, bacteria, and fungi (mycotics)' to bridge the two terms. |
| low | ease-of-use | Line 138 | WDL is introduced as 'a standard workflow language (WDL)' which is acceptable, but 'Docker images' is used without explanation for a non-technical audience. | Briefly gloss 'Docker images' (e.g., 'standardized, reproducible software environments') on first use. |

### `docs/workflows/data_export/microreact_export.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | grammar | Line 9 | Missing comma before a dependent clause and slightly awkward phrasing: 'If an access token is unavailable users can manually upload' needs a comma after 'unavailable'. | Change to 'If an access token is unavailable, users can manually upload the output project file to Microreact.' |
| medium | ease-of-use | Line 9 | Missing 'why/when to use' context and undefined term. The page does not explain what Microreact is or why a public health scientist would export to it; 'Access Token' is capitalized and used without explaining how to obtain one. | Add a sentence describing Microreact (an interactive platform for visualizing genomic data on maps/trees) and briefly note where to get an Access Token (or link to instructions). |
| low | grammar | Line 9 | Compound modifier missing a hyphen: 'phylogenetic trees resulting from other Terra run workflows' reads ambiguously ('Terra run workflows'). | Hyphenate as 'Terra-run workflows' or reword to 'workflows run on Terra'. |
| low | ease-of-use | Lines 16-21 ('Updating Projects' tip) | The example URL 'https://microreact.org/project/&lt;project_url&gt;run' is confusing: it appears to concatenate the placeholder with the word 'run' with no separator, making it unclear what the real project_url looks like. | Clarify the URL pattern and what portion is the project_url (e.g., show a realistic example and highlight exactly which segment to copy). |
| low | ease-of-use | Lines 16-22 | Inconsistent terminology for the same identifier: the text refers to 'project_url', 'project id', and 'The project id' interchangeably, which can confuse a non-technical reader about whether these are the same value. | Use one consistent term, clarifying that the project URL and project ID refer to the same identifier needed for 'project_url'. |

### `docs/workflows/data_export/transfer_column_content.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Lines 13-17 ('Ensure Proper Permissions' tip) | Undefined jargon for non-technical readers: 'Terra pet-service account', 'Terra PROXY account', and 'Storage Object Admin Google privileges' are used without explanation, and the footnote defining GS URI (line 32) is referenced via '[^1]' but the marker is not visible near where GS URIs would be discussed. | Briefly define 'pet-service account' / 'proxy account' or link to a glossary, and ensure the GS URI footnote marker appears at the point where GS URIs are first mentioned. |
| low | grammar | Line 17 | Missing terminal period at the end of the second bullet: '...(Storage Object Admin Google privileges)'. | Add a period for consistency with the first bullet. |

### `docs/workflows/data_import/assembly_fetch.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 12 (Assembly_Fetch_PHB intro, item 2) | The two run options end with inconsistent punctuation (item 1 ends with a comma, item 2 with a period) and item 2 is ambiguous: 'You can provide an assembly' is unclear after item 1 already mentioned 'an accession for the specific assembly'. | Make both list items end with periods and clarify item 2, e.g. 'You can provide your own assembly (FASTA file), and Assembly_Fetch will first use ReferenceSeeker to find the closest reference genome in RefSeq, then download that reference.' |
| medium | ease-of-use | Line 12 / Line 20 ('RefSeq') | 'RefSeq' is used without expansion. For a non-technical audience this NCBI database name is undefined jargon. | Expand on first use, e.g. 'RefSeq (NCBI's Reference Sequence database)'. |
| low | ease-of-use | Line 20 (Inputs) | A single dense run-on sentence describes the required input and the two alternative inputs, which is hard for a non-technical reader to parse. | Break into a short bulleted list: required 'samplename', plus one of 'ncbi_accession' OR 'assembly_fasta'. |

### `docs/workflows/data_import/basespace_fetch.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 62 (Step 2, item 2) | Sentence ends without punctuation: 'the following lines are the commands to be copy/pasted into the terminal' has no closing period. | Add a period at the end of the sentence. |
| low | grammar | Line 149 (Nested Samplenames warning) | Inconsistent product-name capitalization: 'Basespace' is used here while the rest of the document uses 'BaseSpace'. | Change 'Basespace dataset' to 'BaseSpace dataset'. |
| low | ease-of-use | Line 17 (Retrieving BaseSpace Access Credentials) | Comma splice and muddled logic: 'This can be set up in Terra, however it will work in any command-line environment ...' | Split into two sentences: 'This can be set up in Terra. However, it will also work in any command-line environment that has internet access to install and run the BaseSpace command-line tool (bs).' |
| low | ease-of-use | Line 82 (bash code comment) | The comment 'replace the path below with the appropriate path returned by the find command above' references a find command that does not appear earlier in the block, which will confuse users. | Remove the reference to the non-existent find command, or add the find command to the block. |

### `docs/workflows/data_import/create_terra_table.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | accessibility | Lines 88-92 (Outputs table) | The table mapping columns to true/false combinations uses only check/cross emoji to convey meaning, with no textual key. Screen-reader and colorblind users may not get the meaning. | Add a one-line legend, e.g. 'check mark = column is created; cross = column is not created.' |
| low | ease-of-use | Line 81 (Outputs, first bullet) | The phrase 'under the `new_table_name`_id column' references the input variable new_table_name, which is never introduced earlier on the page. | Briefly note that new_table_name is the input where the user specifies the desired table name, or link to its row in the Inputs table. |

### `docs/workflows/data_import/fetch_srr_accession.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Lines 9-11 (intro) | The intro assumes the reader understands the relationship between BioSample IDs, SRA Experiment IDs and SRR accessions, and gives no 'why' for converting one to the other. | Add one sentence on the purpose, e.g. 'An SRR accession is needed to download the actual sequencing reads (for example, with the SRA_Fetch workflow). This workflow finds the SRR when you only have a BioSample or Experiment ID.' |
| low | accessibility | Lines 28-33 (techdetails table) | The techdetails admonition and table appear under-indented relative to the enclosing '??? task' block (table starts at column 0 while the admonition is at 4 spaces), which may break rendering inside the toggle. | Verify the techdetails admonition and table are indented consistently inside the task toggle so the table renders correctly. |

### `docs/workflows/data_import/ont_barcode_concatenation.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Title / Line 7-9 (intro) | 'ONT' is used in the title and throughout but is never expanded to Oxford Nanopore Technologies. A non-technical reader may not know what ONT means. | Expand on first use, e.g. 'This workflow concatenates Oxford Nanopore Technologies (ONT) sequencing reads ...'. |
| low | accessibility | Line 103 ('clicking on the two squares') | A UI control is identified only by its shape ('the two squares next to the file path'), which is not helpful for screen-reader users or those who interpret the icon differently. | Add the control's function, e.g. 'clicking the copy-path icon (two overlapping squares) next to the file path'. |
| low | ease-of-use | Lines 9, 39 (terminology) | The page alternates between 'folder', 'directory', and 'subdirectory' for the same concept, which can confuse a non-technical reader. | Standardize on one term (e.g., 'folder') or state explicitly that 'folder' and 'directory' mean the same thing. |
| low | ease-of-use | Lines 21-22 (directory-structure diagram) | The input tree shows a stray 'output_bucket_path/' line at the top with no contents and no explanation, even though the surrounding text discusses only input_bucket_path. This makes the diagram ambiguous. | Remove the stray 'output_bucket_path/' line from the input tree (it is shown separately later) or annotate that it is a separate, empty folder. |

### `docs/workflows/data_import/sra_fetch.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 35 (Outputs) | 'Given the lack of usefulness of SRA Lite formatted FASTQ files ...' assumes the reader understands why artificial (Q30/Q3) quality scores are problematic. | Add a brief clause on the impact, e.g. 'Because SRA Lite files have artificial quality scores, they are less useful for quality-dependent analyses; we therefore try to avoid them ...'. |
| low | grammar | Line 9 (intro) | Run-on with wrong preposition: 'It requires an SRA run accession then populates the associated read files to a Terra data table.' | Reword, e.g. 'It requires an SRA run accession, and then populates a Terra data table with the associated read files.' |
| low | grammar | Line 15 (Inputs) | Missing word: 'an SRA run accession beginning "SRR"' should be 'beginning with "SRR"' (and likewise for ERR and DRR). | Change to 'beginning with "SRR" ... beginning with "ERR" ... beginning with "DRR".' |
| low | accessibility | Lines 19-24 (accession naming-scheme table) | The table has no header row; it relies on the reader inferring that the first column is object type and the second is the accession prefix, which hampers screen-reader navigation. | Add header cells, e.g. 'Object Type \| Accession Prefix'. |

### `docs/workflows/genomic_characterization/freyja.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | typo | Line 34 (Freyja_Dashboard_PHB bullet) | Stray trailing '=' character after the link: '[**Freyja_Dashboard_PHB**](freyja.md#freyja_dashboard)='. | Remove the trailing '=' so the line ends with the link. |
| medium | ease-of-use | Lines 9-12, 50 (Overview / Data Quality) | 'Demixing' and 'deconvolution' are core to the tool but are never defined in plain language for a non-technical reader. | Add a brief plain-language definition of 'demixing' on first use (separating a mixed sample into the proportion of each viral lineage present). |
| low | typo | Line 158 (freyja_demixed format table) | Mismatched quote in the example value: "[('Delta', 0.65), ('Other', 0.25), ('Alpha', 0.1')]" has a stray apostrophe after 0.1. | Change '0.1'' to '0.1' to fix the mismatched quote. |
| low | grammar | Line 60 (Freyja_FASTQ_PHB Inputs intro) | Wrong article: 'compatible with the multiple input data types'. | Change to 'compatible with multiple input data types'. |
| low | grammar | Line 165 (Outputs, lineage array bullet) | Inconsistent field name: the bullet says 'The `lineage` array' but the table and surrounding text use 'lineages'. | Change '`lineage` array' to '`lineages` array'. |
| low | accessibility | Lines 156-162 (freyja_demixed format table) | The transposed key/value table has a blank first header cell and is shown without explanation of its unusual layout, hampering screen-reader interpretation. | Add a sentence before the table noting it lists row labels (summarized, lineages, abundances, resid, coverage) in the first column with the corresponding values for one sample in the second. |
| low | ease-of-use | Lines 50-52 ('Freyja, Sequencing Platforms and Data Quality') | Two long, jargon-heavy paragraphs ('deconvolution process', 'Q-scores', '100X coverage', 'cocirculating lineages') with a near-duplicate sentence about sequencing depth increasing as quality decreases. | Break into bullet points, briefly define Q-score and coverage, and remove the near-duplicate sentence. |
| low | ease-of-use | Lines 311-322 vs 343-354 (two organism lists) | The barcode-supported organism list and the freyja_pathogen allowed-options list use different names and partly different organisms with no explanation of the difference, confusing users choosing a value. | Add a note explaining that one list is for the freyja_pathogen flag and the other is for available barcode files, and why the names/sets differ. |

### `docs/workflows/genomic_characterization/pangolin_update.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | grammar | Line 9 (intro) | Garbled sentence: 'The Pangolin_Update workflow re-runs Pangolin updating prior lineage calls from one docker image to meet the lineage calls specified in an alternative docker image.' | Reword, e.g. 'The Pangolin_Update workflow re-runs Pangolin to update previous lineage calls, replacing the assignments made by an older Pangolin docker image with those of a newer one.' |
| medium | ease-of-use | Line 9 ('docker image') | 'docker image' is the central concept but is never explained; a non-technical reader may not know what it is or why lineage calls depend on it. | Add a plain-language note, e.g. 'A docker image bundles a specific version of the Pangolin software and its lineage database; using a newer image updates your lineage assignments to the latest nomenclature.' |
| low | ease-of-use | Lines 7-9 (page intro) | The page does not clearly state when to use this workflow beyond one embedded clause. | Add a sentence such as 'Use this workflow when SARS-CoV-2 lineage nomenclature has been updated and you want to refresh older results without re-running the full characterization workflow.' |

### `docs/workflows/genomic_characterization/theiacov.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | accessibility | Lines 458-473 (Organism-specific characterization table) | The table uses three glyph-only symbols (check, plus, cross) whose legend appears only after the table, so screen-reader users navigating cell-by-cell encounter the symbols before the key. | Place the symbol legend immediately before the table or as a caption so the meaning is available first. |
| medium | ease-of-use | Line 13 and throughout (acronyms ONT, VADR, IRMA, GenoFLU) | Several acronyms are used without expansion on first use (ONT, VADR, IRMA, GenoFLU), which is a problem for the explicitly non-technical audience. | Expand acronyms on first use (e.g., ONT = Oxford Nanopore Technologies) or add a brief note that the items in the characterization table are tool names. |
| low | ease-of-use | Line 86 | Dangling forward reference: 'We've provided the following information ... in the form of input JSONs.' appears at the end of the Supported Organisms section, but the JSON links are above in the Key Resources sidebar (lines 29-51), not below. | Point to the Key Resources panel explicitly, e.g. 'See the input JSONs in the Key Resources panel above.' |
| low | ease-of-use | Lines 218-219 (West Nile Virus table) | The value 'NA' is used in tables to mean 'not applicable/unsupported' but this convention is never stated, so a non-technical reader may misread it. | Add a brief general note that 'NA' in these tables means the parameter is not used for that organism. |

### `docs/workflows/genomic_characterization/theiaeuk.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 25 (Before running TheiaEuk warning) | Redundant/awkward wording: 'the read support reported is, at the moment, non-reliable. Future improvements will include improvements on this module.' ('non-reliable' is non-idiomatic; 'improvements ... improvements' repeats). | Reword, e.g. 'the reported read support is currently unreliable. We plan to improve this module in the future.' |
| low | grammar | Line 57 (Workflow Tasks intro) | Subject-verb agreement: 'TheiaEuk workflow subsequently launch default genome characterization modules' (singular subject, plural verb). | Change to 'The TheiaEuk workflows subsequently launch ...' or 'workflow subsequently launches ...'. |
| low | ease-of-use | Line 62 (Core tasks note) | Awkward, hard-to-parse phrasing duplicated from the intro: 'tasks that are performed regardless of and specific for the input data type'. | Reword, e.g. 'These core tasks run for every sample. Some are common to all data types, while others are specific to the input data type (Illumina or ONT).' |
| low | ease-of-use | Line 11 ('in silico') | '_in silico_' is used without explanation for a non-technical reader (AMR is expanded, which is good). | Briefly gloss 'in silico' as 'computational / by software' on first use. |

### `docs/workflows/genomic_characterization/theiameta.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Lines 9, 13-14 (intro) | Dense with undefined terms for a non-technical audience: 'contigs', 'de novo', 'binning', 'taxon'. The binning sentence assumes the reader already knows what contigs and taxa are. | Add brief glosses for 'contig' (a contiguous assembled sequence) and 'de novo' (assembled from scratch, without a reference) on first use. |
| low | grammar | Line 9 (intro) | Awkward phrasing and inconsistent hyphenation: 'either using a reference-genome or not' ('reference-genome' is hyphenated here but 'reference genome' elsewhere). | Reword to 'with or without a reference genome' and remove the hyphen. |
| low | grammar | Lines 64-65 and 146-147 (samtools task) | Trailing space in the heading 'samtools: SAM File Conversion ' and repeated verb: 'converts the output SAM file from minimap2 and converts it to a BAM file' (double 'converts', also at line 147). | Remove the trailing space and reword to 'takes the output SAM file from minimap2 and converts it to a BAM file'. |
| low | ease-of-use | Line 44 (metaspades task description) | Dense, jargon-heavy paragraph ('de Bruijn graph', 'graph simplification procedures', 'highly polymorphic diploid genomes') with no plain-language summary. | Add a one-sentence plain-language summary (metaSPAdes assembles short reads into longer genome fragments) before the technical detail. |

### `docs/workflows/genomic_characterization/theiaprok.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 152 / Line 164 / Line 229 (acronyms AMR, MLST, XDR) | Multiple acronyms used without expansion for a non-technical audience: AMR ('AMR Characterization'), MLST (Sequence Type section), and XDR (line 229). | Expand AMR (antimicrobial resistance), MLST (multi-locus sequence typing), and XDR (extensively drug-resistant) on first use. |
| low | ease-of-use | Line 84 (Core Tasks intro) | Awkward, hard-to-parse phrasing: 'tasks that are performed regardless of and specific for the input data type'. | Reword, e.g. 'Some core tasks run for all data types, while others are specific to the input data type.' |
| low | ease-of-use | Line 229 (Shigella XDR prediction tip) | 'Shigella XDR prediction' uses the undefined acronym XDR and only vaguely points to 'the documentation section above for ResFinder'. | Expand XDR (extensively drug-resistant) and briefly state what the ResFinder section explains. |

### `docs/workflows/genomic_characterization/theiaviral.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 99 (TheiaViral_Panel input section) | Inconsistent capitalization/styling of 'skani' vs 'Skani' across the page (e.g. line 99 'best retrieved skani reference', line 408 'custom Skani database', line 424 'reference selection software, Skani'). Inconsistent terminology for the same tool can confuse non-technical readers. | Standardize on one form (e.g. 'Skani') throughout the document. |
| medium | ease-of-use | Line 370 (Taxa-Specific Tasks intro) | The phrase 'using the taxon ID integer because these can be simpler' is vague for a non-technical reader and does not explain why an integer ID is simpler than a name. | Clarify, e.g. 'We recommend using the numeric taxon ID (e.g. 2697049) because it avoids spelling and capitalization mismatches that can occur with full names.' |
| low | typo | Line 126 (heading) and Line 127 (`extract_unclassifed` ... set to `false`) | Within the `extract_unclassified` toggle block, the parameter name is spelled `extract_unclassifed` (missing the second 'i') in the body text on line 127, which is inconsistent with the heading and other references. | Correct the body text to `extract_unclassified` to match the parameter name used elsewhere. |
| low | grammar | Line 113 (`output_taxon_table` toggle) | Subject-verb / number agreement: 'specifies which taxon are output to what taxon table' uses singular 'taxon' with plural 'are'. | Rephrase to 'specifies which taxa are output to which taxon table' (and consider 'which' instead of 'what'). |
| low | grammar | Line 124 (`kraken_db` toggle) | Missing comma before a dependent clause: 'When making changes to this parameter keep in mind the relationship...' is missing a comma after the introductory clause. | Change to 'When making changes to this parameter, keep in mind the relationship between these two inputs.' |
| low | grammar | Line 451 ('How is consensus assembly quality evaluated?') | Missing article: 'if the reference genome is sufficient quality' should read 'of sufficient quality' (this wording also appears on line 418). | Change to 'if the reference genome is of sufficient quality'. |
| low | grammar | Line 481 (`metaviralspades_status` bullet) | Awkward/incomplete phrase: 'extract viral contigs in the _de novo_' on line 479 reads 'submitting the _de novo_ assembly'; on line 481 'fallback _de novo_ assemblers are implemented for both TheiaViral workflows' is fine, but line 478 'BLAST whatever contigs do exist in the _de novo_ to determine' omits the noun 'assembly' after '_de novo_'. | Add the noun: 'BLAST whatever contigs do exist in the _de novo_ assembly to determine...'. |
| low | accessibility | Line 19 (tip block, Step 2) and first mention of ANI | The acronym ANI is introduced as 'average nucleotide identity (ANI)' which is good, but 'RefSeq' (line 51) and 'gVCF'-style terms aside, 'RefSeq _assembly_ accessions' uses 'RefSeq' without expansion on first use for a non-technical audience. | Briefly gloss 'RefSeq' on first use (e.g. 'RefSeq, NCBI's curated Reference Sequence collection') or link to a definition. |
| low | ease-of-use | Line 104 (`taxon_ids` toggle) | Redundant/repetitive instruction: the passage states twice in nearly identical wording that the taxon IDs must be present in the Kraken2 database ('the taxon IDs _must_ be present in the Kraken2 database' then 'Keep in mind that these IDs must be available in the passed Kraken DB'). | Remove the duplicate sentence to tighten the paragraph. |
| low | ease-of-use | Line 104 and elsewhere (Kraken2 / Kraken DB) | Inconsistent terminology: 'Kraken2 database', 'Kraken DB', and 'Kraken database' are used to refer to the same thing, which may confuse non-technical readers. | Use a single consistent term, e.g. 'Kraken2 database', throughout. |

### `docs/workflows/genomic_characterization/vadr_update.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | accessibility | Line 9 ('VADR_Update_PHB' section, first paragraph) | The acronym VADR is never expanded on first use. A non-technical reader will not know what VADR stands for or what it does. | Expand on first use, e.g. 'VADR (Viral Annotation DefineR), a tool that validates and annotates viral genome sequences', and add a one-line statement of when/why a user would run this workflow. |
| medium | ease-of-use | Line 7-9 (VADR_Update_PHB section) | Missing 'why/when to use this' guidance. The page explains what the workflow runs (VADR with separately provided models) but not the purpose of validating/annotating sequences or when a public health scientist would choose this standalone workflow. | Add a sentence describing the purpose (e.g. checking that consensus genomes meet submission/annotation quality standards) and when to use VADR_Update versus running VADR inside another workflow. |
| low | accessibility | Line 11-21 (organism/model table) | The recommended-models table is presented with no textual lead-in explaining how to read it or how to use the four columns (vadr_model_file, vadr_opts, max_length). Screen-reader users and non-technical readers get a dense table of paths and flags without context. | Add a sentence before the table explaining that each row gives the model file, recommended command-line options, and maximum sequence length to set for that organism's input parameters. |
| low | ease-of-use | Line 9 | Jargon 'slimmed-down docker image' and 'requires models to be provided separately' may not be clear to a non-technical reader, and it is not stated where to obtain the models (the table that follows provides them but this connection is implicit). | Add a clause pointing the reader to the table below for the model files to provide, e.g. 'you must supply a model file separately (see the table below for recommended models per organism).' |

### `docs/workflows/phylogenetic_construction/augur.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 39 (geographical information bullet, 'divisions') | Inconsistent/confusing level naming: the resolution levels are described as 'Lowest-level resolution' for region (continents) and 'Highest-level resolution' for location (cities). This is counterintuitive—continents are the broadest (lowest detail) and cities the finest (highest detail)—but a non-technical reader may read 'lowest-level' as 'least important' or 'smallest area'. | Clarify wording, e.g. 'region — broadest resolution (least geographic detail), used for continents' and 'location — finest resolution (most geographic detail)'. |
| low | grammar | Line 34 (How to prepare metadata bullet) | Number agreement: 'You can specify unknown dates or month' mixes plural 'dates' with singular 'month'. | Change to 'unknown dates or months'. |
| low | grammar | Line 34 (How to prepare metadata bullet) | Awkward phrasing: 'by replacing the respective values by `XX`' uses 'by' where 'with' is expected. | Change to 'by replacing the respective values with `XX`'. |
| low | grammar | Lines 78 and 85 (metadata logic paragraph and caption) | Missing article: 'metadata that is present in input metadata file' should be 'in the input metadata file'. This appears on both line 78 and the duplicated caption text on line 85. | Insert 'the': 'present in the input metadata file'. |
| low | grammar | Line 106 ('Running Augur_PHB on custom organisms') | Logic/wording error: 'For non-default organisms (listed above)' is contradictory because the organisms listed above are the DEFAULT organisms; non-default organisms are those NOT listed. | Change to 'For organisms other than those listed above' or 'For non-default organisms (i.e. not listed above)'. |
| low | grammar | Line 119 | Missing terminal punctuation: the sentence ending '...detailed documentation on Augur](...)' has no period. | Add a period at the end of the sentence. |
| low | grammar | Line 106 | Missing terminal punctuation: 'several optional inputs are required to guarantee workflow functionality' ends without a period. | Add a period at the end of the sentence. |
| low | accessibility | Lines 9-10 (Augur Workflow Overview) and 80-86 (How Metadata Shapes the Augur Output) | Augur and Auspice are core tools but their relationship/roles are introduced somewhat far into the text; 'JSON' (line 82, 159) is used without expansion for a non-technical audience. | Briefly gloss 'JSON' on first use as a file format (e.g. 'an Auspice JSON file, a structured data file Auspice reads'). |
| low | ease-of-use | Line 39 (geographical information bullet) | The level is named 'divisions' (plural) in the bullet but Augur fields are typically singular ('division'); inconsistent with 'region', 'country', 'location' which are singular, which could confuse users entering the field name. | Verify and use the exact field name consistently (singular 'division' if that is the field name). |

### `docs/workflows/phylogenetic_construction/clair3_variants.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 9 (Clair3_Variants_ONT first paragraph) | Dense, jargon-heavy single paragraph for a non-technical audience: 'deep learning-based variant calling', 'sorted and indexed BAM files', 'gVCF format' are introduced without definition, and the paragraph is long. | Break into shorter sentences/bullets and briefly gloss key terms (BAM, gVCF) or add a one-line plain-language summary of what the workflow produces and why. |
| low | accessibility | Line 9 and line 89 (gVCF / structural variants) | Acronyms used without expansion on first use: 'gVCF' (line 9) and 'SNVs' (line 87, distinct from 'SNPs' defined on line 9). The body also lists 'Structural variants' (line 89) without explanation. | Expand 'gVCF' (genome VCF) and 'SNVs' (single nucleotide variants) on first use, and note that SNVs/SNPs are used interchangeably here to avoid confusion. |
| low | ease-of-use | Lines 32-34 (Supported Clair3 Models table) and line 9 | Inconsistent terminology between table and text: line 9 says it detects 'single nucleotide polymorphisms (SNPs)' while line 87 says 'Single nucleotide variants (SNVs)'. Using two terms for the same concept can confuse non-technical readers. | Standardize on one term (SNPs or SNVs) or explicitly state they are equivalent. |
| low | ease-of-use | Lines 41-46 (Note on Haploid Settings) | Assumes the reader knows what 'haploid', 'phasing', and 'homozygous variants (1/1)' mean. A non-technical public health scientist may not, and the note does not explain when these defaults should be changed. | Add a brief plain-language note on what haploid means (single copy of the genome, typical for bacteria/viruses) and when a user might need to adjust these defaults. |

### `docs/workflows/phylogenetic_construction/core_gene_snp.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 14 (Core_Gene_SNP_PHB intro) | The intro states what tools are used (PIRATE, IQ-TREE) but does not explain in plain language what 'core genome', 'pangenome', or 'SNP' analysis is, or when/why a user would choose this workflow. Audience is non-technical. | Add a sentence defining core genome vs pangenome and describing the use case (e.g. comparing closely related bacterial isolates for outbreak analysis). |
| low | accessibility | Line 14 (first use of SNP and GFF3) | 'SNP' appears in the title/intro without expansion on this page, and the workflow diagram alt text references 'GFF3 genome annotation files' (a file-format term) that is not explained in the body for non-technical readers. | Expand SNP (single nucleotide polymorphism) on first use and add a brief note that GFF3 files are genome annotation files (e.g. produced by an annotation workflow). |
| low | ease-of-use | Lines 33-48 (Core Tree vs Pangenome Tree generation blocks) | The 'Core Tree Generation' and 'Pangenome Tree Generation' blocks include the identical set of task includes (snp_sites, IQ-TREE, snp_dists, reorder_matrix), which may appear duplicated/confusing without a note explaining they are the same tasks applied to different alignments. | Add a brief note clarifying that the same set of tasks runs on the core-genome alignment and (optionally) the pangenome alignment, controlled by `core_tree` and `pan_tree`. |

### `docs/workflows/phylogenetic_construction/find_shared_variants.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 9 (Find_Shared_Variants_PHB intro) | Missing 'why/when to use this' guidance for a non-technical reader: the description explains the mechanics (concatenating and reshaping variant results) but not the purpose, e.g. identifying mutations shared across an outbreak or sample group. | Add a sentence describing the use case, e.g. 'This is useful for identifying mutations common to a group of related samples, such as isolates from a suspected outbreak.' |

### `docs/workflows/phylogenetic_construction/ksnp3.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | typo | Line 16 (video link) | Inconsistent capitalization of the tool name: link text reads 'Using KSNP3' (all caps KSNP) while the rest of the page and the actual tool name use 'kSNP3'. | Change 'KSNP3' to 'kSNP3' in the link text for consistency. |
| low | grammar | Line 14 (summarize_data description) | Missing comma in a list and missing comma before 'else': 'such as AMR genes, plasmid types etc.' should have a comma before 'etc.', and 'for visualization in Phandango), else it can be viewed in Excel' uses 'else' as a conjunction awkwardly. | Use 'such as AMR genes, plasmid types, etc.' and rephrase to '...for visualization in Phandango; otherwise it can be viewed in Excel.' |
| low | accessibility | Line 14 (summarize_data description) | Acronym 'AMR' is used without expansion on first use. | Expand on first use: 'AMR (antimicrobial resistance) genes'. |
| low | accessibility | Line 16 | The emoji 'TV' (📺) prefixes the video link; emoji can be read aloud unhelpfully by screen readers and conveys meaning by icon alone. The link text itself is descriptive, so the emoji is decorative. | Consider prefacing with the word 'Video:' rather than relying on the emoji, or ensure the descriptive link text stands alone. |

### `docs/workflows/phylogenetic_construction/ksnp4.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 14 (summarize_data description) | Same as ksnp3: 'such as AMR genes, plasmid types etc.' lacks a comma before 'etc.', and 'else it can be viewed in Excel' uses 'else' awkwardly as a conjunction. | Use 'such as AMR genes, plasmid types, etc.' and rephrase to '...for visualization in Phandango; otherwise it can be viewed in Excel.' |
| low | accessibility | Line 14 (summarize_data description) | Acronym 'AMR' is used without expansion on first use. | Expand on first use: 'AMR (antimicrobial resistance) genes'. |
| low | accessibility | Line 16 | The emoji 'TV' (📺) prefixes the video link; it conveys meaning by icon alone and may be read unhelpfully by screen readers (the descriptive link text is fine on its own). | Consider prefacing with 'Video:' instead of relying on the emoji. |

### `docs/workflows/phylogenetic_construction/lyve_set.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 35 (SNP discovery step) | Double space before the parenthetical: 'from the Lyve-SET pipeline  (`smalt` and `varscan`)' contains two spaces. | Remove the extra space. |
| low | grammar | Line 57 (Outputs section) | Missing word: 'These files can be accessed viewing the execution directory' is missing 'by'. | Change to 'These files can be accessed by viewing the execution directory for the run.' |
| low | accessibility | Line 9 (Lyve_SET_PHB intro) | Acronym 'hqSNPs' is introduced as 'high quality single nucleotide polymorphisms (hqSNPs)' which is good, but 'RAxML' is named on line 9 with no expansion or link, and 'WDL' is used without expansion (line 9, 28). | Briefly note that RAxML is a maximum-likelihood phylogenetics tool (and/or link it), and expand 'WDL' (Workflow Description Language) on first use. |
| low | ease-of-use | Line 31 (Read processing step) | Jargon 'CG-Pipeline ("CGP")' and 'BayesHammer' are named as read-cleaning options without any explanation of what they do or how to choose between them, which a non-technical reader cannot evaluate. | Add a one-line description of each option or guidance on when to choose one over the other (or default). |
| low | ease-of-use | Line 33 (Reference procurement step) | Jargon 'mask phages or cliffs' is partially explained (cliffs are defined) but 'phages' and the rationale 'remove low quality SNPs' assume background knowledge; non-technical readers may not know what masking a phage region accomplishes. | Briefly note why one might mask phage regions (e.g. they are mobile/variable regions that can introduce misleading SNPs). |

### `docs/workflows/phylogenetic_construction/mashtree_fasta.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 11 (MashTree_FASTA_PHB intro paragraph) | Dense, technical explanation aimed above the stated audience: 'kmers', 'integer value with hashing and Bloom filters', 'sketch', and 'neighbor-joining algorithm' are all introduced in one paragraph without plain-language definitions. | Break into shorter sentences and add a brief plain-language gloss for 'k-mer' and 'neighbor-joining', or add a one-line summary of what the workflow does and when to use it (fast approximate trees from assemblies). |
| low | accessibility | Line 11 (kmers) | Term 'kmers' is used repeatedly without expansion/definition; standard styling is 'k-mers' and the concept is undefined for non-technical readers. 'AMR' (line 13) is also unexpanded. | Define 'k-mer' on first use (short subsequences of length k) and expand 'AMR (antimicrobial resistance)'. |
| low | ease-of-use | Line 9 / Line 11 | Missing 'when to use this' guidance: the page explains how Mash distances are computed but not why a user would choose MashTree (e.g. very fast, approximate clustering of many genomes) over SNP-based methods. | Add a sentence on the use case and tradeoff (speed vs resolution) relative to SNP-based phylogenetic workflows. |

### `docs/workflows/phylogenetic_construction/snippy_streamline.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 22 (Reference Genome Options, option 3) | Long, dense sentence describing automatic reference selection ('automatically selected using the `centroid` task and `reference_seeker` task to find a close reference genome to your dataset by providing data in the `assembly_fasta` input variable and leaving the ... fields blank'). It is hard to parse for a non-technical reader. | Break into a shorter sentence or sub-bullets separating the action (auto-select) from the required input (provide `assembly_fasta`) and the condition (leave the other two fields blank). |
| low | grammar | Line 41 (Inputs intro) | Missing terminal punctuation: the sentence ending '...the closest reference genome to your dataset (via `assembly_fasta`)' has no period. | Add a period at the end. |
| low | grammar | Line 32 (Phylogenetic Tree Construction Options, option 3) | Missing noun: 'choosing the nucleotide substitution (_by altering `iqtree2_model`...)' omits 'model' after 'substitution'. | Change to 'choosing the nucleotide substitution model'. |
| low | accessibility | Line 66 (MTBC) and line 64 heading | Acronyms 'MTBC' and 'MTB' are used (line 66) without expansion in the body; 'MLST' (line 50) is also unexpanded. The heading expands 'Mycobacterium tuberculosis complex' but the body acronym 'MTBC' is not tied to it explicitly for all readers. | Define MTBC on first use ('Mycobacterium tuberculosis complex (MTBC)') and expand MLST (multilocus sequence typing). |
| low | ease-of-use | Lines 110-111 (Gubbins Nucleotide Substitution Model tip) and throughout | Inconsistent capitalization of tool name 'gubbins' vs 'Gubbins' (e.g. line 33 'detected by gubbins', line 37 'incompatible with Gubbins', line 111 'used by gubbins'). | Standardize on 'Gubbins' throughout. |
| low | ease-of-use | Lines 32, 60 (IQ-TREE naming) | Inconsistent tool naming: 'IQ-TREE' (line 32), 'IQTree's Model Finder' (line 62), and 'IQTREE' in references vary across the page; 'ModelFinder' vs 'Model Finder' also varies. | Standardize on 'IQ-TREE 2' and 'ModelFinder' consistently. |

### `docs/workflows/phylogenetic_construction/snippy_streamline_fasta.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 32 (Phylogenetic Tree Construction Options, option 3) | Missing noun: 'choosing the nucleotide substitution (_by altering `iqtree2_model`...)' omits 'model' after 'substitution'. | Change to 'choosing the nucleotide substitution model'. |
| low | accessibility | Line 64 (heading) and line 64-73 body (MTBC/MTB/MLST) | Acronyms 'MTBC'/'MTB' (line 64) and 'MLST' (line 48) are used without expansion in the body text. | Define MTBC on first use ('Mycobacterium tuberculosis complex (MTBC)') and expand MLST (multilocus sequence typing). |
| low | ease-of-use | Lines 33, 37, 109 (Gubbins capitalization) | Inconsistent capitalization of tool name 'gubbins' vs 'Gubbins' (e.g. line 33 'detected by gubbins', line 37 'incompatible with Gubbins', line 109 'used by gubbins'). | Standardize on 'Gubbins' throughout. |
| low | ease-of-use | Line 24 (Reference Genome Options, option 3) | The automatic reference selection option does not state what input the user must provide (unlike snippy_streamline.md which mentions `assembly_fasta`). For this FASTA workflow the assemblies are the main input, but the omission could leave a non-technical reader unsure of what enables auto-selection. | Add a clause clarifying that the input assemblies (FASTA) are used to find the closest reference when the reference fields are left blank. |

### `docs/workflows/phylogenetic_construction/snippy_tree.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 9, Snippy_Tree_PHB section | The acronym SNP is used ("pairwise SNP-distance matrix") without being expanded on first use. A non-technical reader may not know SNP means single-nucleotide polymorphism. | Expand on first use, e.g. "pairwise single-nucleotide polymorphism (SNP) distance matrix". |
| low | grammar | Line 35, 'Default settings' toggle | Sentence ends without a period: "...suitable for generating phylogenies for most bacteria". | Add a closing period. |
| low | grammar | Line 43, '_Mycobacterium tuberculosis_ complex' toggle | Sentence "Phylogenies of MTBC are typically constructed" lacks ending punctuation before the bulleted list that completes it; the colon/lead-in is missing. | End with a colon: "Phylogenies of MTBC are typically constructed:" |
| low | accessibility | Line 43, '_Mycobacterium tuberculosis_ complex' toggle | MTBC acronym introduced ("Phylogenies of MTBC are typically constructed") with no expansion; the heading uses the full name but the body abbreviates it, which a screen-reader/non-technical user may not connect. | Write "Phylogenies of the Mycobacterium tuberculosis complex (MTBC) are typically constructed..." |
| low | ease-of-use | Line 14, Snippy_Tree_PHB section (bullet list) | "with a bed file" and "TB" are used without explanation. "bed file" is jargon for a genomic region format, and TB (tuberculosis) is unexpanded. | Briefly note that a BED file specifies genomic regions to mask, and expand TB to tuberculosis on first use. |
| low | ease-of-use | Line 16, bullet list | "Decide which nucleotide substitution model to use" gives no guidance on what this means or when it matters for a non-technical user. | Add a brief note (or link) explaining that the default automatically selects a model, so most users do not need to change it. |

### `docs/workflows/phylogenetic_construction/snippy_variants.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 19, 'Assessing support for a mutation' bullet | IGV is used as an acronym without expansion ("can be visualized in IGV"). Non-technical readers will not know IGV is the Integrative Genomics Viewer. | Expand on first use: "visualized in IGV (Integrative Genomics Viewer)". |
| low | grammar | Line 20, 'Example Use Cases' | "may be an error of sequencing" is awkward phrasing. | Reword to "may be a sequencing error." |
| low | ease-of-use | Line 19, 'Assessing support for a mutation' bullet | "see Theiagen Office Hours recordings" provides no link or guidance on where to find these recordings. | Add a link to the Office Hours recordings or state where they are located. |
| low | ease-of-use | Line 9, Snippy_Variants_PHB section | Acronyms SNPs, MNPs, INDELs are introduced together; while SNP/MNP/INDEL are spelled out once, the dense parenthetical "(SNPs), multi-nucleotide polymorphisms (MNPs), and insertions/deletions (INDELs)" in one long sentence may be hard for non-technical readers. | Consider a short bulleted list of the variant types to improve readability. |

### `docs/workflows/phylogenetic_placement/nextclade_batch.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 11, Nextclade_Batch_PHB section | "mutation-annotated tree" is referenced and the term "Augur" introduced without context; a non-technical reader will not know what these are. | Briefly define mutation-annotated tree and note that generating one is an advanced task most users will not need. |
| low | accessibility | Line 11, Nextclade_Batch_PHB section | Link text "the Augur documentation" is descriptive, which is good, but "Contact us" gives no contact method/link. | Add a contact link or email (e.g. support@theiagen.com) after "Contact us". |
| low | ease-of-use | Line 9, Nextclade_Batch_PHB section | The phrase "genotypes batches of samples" and the dense first paragraph mix several concepts (mutation calling, placement, genotyping) in one block. | Add a short "When to use this workflow" note and consider breaking the paragraph into bullets. |

### `docs/workflows/phylogenetic_placement/usher.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 9, Usher_PHB section | Several terms used without explanation for non-technical readers: "maximum parsimony", "protobuf format", and acronyms MAT/UCSC. "Contact us" has no contact method. | Briefly explain maximum parsimony in plain language, note protobuf is a file format, expand UCSC (University of California, Santa Cruz), and add a contact link/email. |
| low | ease-of-use | Line 13, Inputs section | "set-level" and "sample-level" workflow terminology is used without explanation; non-technical Terra users may not know what these mean. | Add a one-line definition or link explaining set-level vs sample-level workflows. |

### `docs/workflows/public_data_sharing/mercury_prep_n_batch.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | accessibility | Lines 16-20, repository support table | The table uses a checkmark symbol (✓) to indicate support with no textual key/legend; screen readers may announce it inconsistently, and there is no text alternative explaining what a checkmark vs blank cell means. | Add a sentence before or after the table stating that a checkmark means the organism/repository combination is supported and a blank cell means it is not. |
| medium | ease-of-use | Line 14, Mercury_Prep_N_Batch_PHB section | Acronyms in the support table (BankIt, BioSample, GenBank, GISAID, SRA) are presented with no explanation of what each repository is or when you would submit to it. | Add a brief one-line description of each repository, or link to each, so users know which target applies to them. |
| low | accessibility | Line 23, 'Mercury expects data tables made with TheiaCoV' admonition | Vague link text "this file" ("See this file for the hard-coded list..."). | Use descriptive link text such as "see the Metadata.py source file". |
| low | accessibility | Line 9, 'Command-line incompatible' admonition | GCP acronym is expanded on line 12 but the repository table and broader doc use acronyms like NCBI, GISAID, SRA, PH4GE that are not all expanded on first use for non-technical readers. | Expand SRA (Sequence Read Archive) and confirm NCBI/GISAID are expanded at first use; PH4GE is expanded on line 12 (good). |
| low | ease-of-use | Lines 64-65, 'Use the sample table' tip | "root entity" and "set table" vs "sample table" are confusing without a quick explanation; the instruction is potentially ambiguous for non-technical Terra users. | Briefly clarify what root entity means and why the set table is the root while terra_table_name should point to the sample table. |

### `docs/workflows/public_data_sharing/terra_2_ena.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 161, Viral Optional Fields, receipt_date row | Missing space after colon: "Format:YYYY-MM-DD." | Add a space: "Format: YYYY-MM-DD." |
| low | accessibility | Line 19, Prerequisites section | Curly/smart quotes used around button names (‘Register Study’, ‘Studies Report’) are inconsistent with the straight quotes used elsewhere and can read oddly with assistive tech. | Use consistent straight quotes or bold/code formatting for UI button names. |
| low | ease-of-use | Line 9, Terra_2_ENA_PHB section | ENA acronym is never expanded (European Nucleotide Archive); appears only as ENA throughout. | Expand on first use: "the European Nucleotide Archive (ENA)". |
| low | ease-of-use | Line 88, Bacterial Metadata Mandatory Fields warning; line 124 Viral Metadata | INSDC acronym used in the missing-values note without expansion (International Nucleotide Sequence Database Collaboration). | Expand INSDC on first use. |
| low | ease-of-use | Line 139, Viral Mandatory Fields, host_sex row | Description says "Gender or sex of the host" which conflates gender and sex; for a metadata field labeled host sex this is imprecise. | Use "Sex of the host" to match the field name and avoid conflating the two terms. |

### `docs/workflows/public_data_sharing/terra_2_gisaid.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 25, 'GISAID Credentials' warning | Inconsistent naming/formatting: "authentication_file" is plain text while "gisaid_credentials" is plain text too, but neither is in code format despite being input variable names referenced elsewhere as code. | Format both as inline code (`authentication_file`, `gisaid_credentials`) for consistency. |
| low | ease-of-use | Line 9, Terra_2_GISAID_PHB section | GISAID acronym is never expanded; non-technical readers may not know it is the Global Initiative on Sharing All Influenza Data. | Expand GISAID on first use or link to its site with a brief description. |
| low | ease-of-use | Line 18, Inputs section | "frameshift" is unexplained jargon for a non-technical reader, and the three option descriptions are quoted verbatim from a web UI without plain-language context on which to choose. | Add a brief explanation of what a frameshift is and a note recommending the default (catch_novel) for most users. |

### `docs/workflows/public_data_sharing/terra_2_ncbi.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | grammar | Line 96, 'What's the difference...' toggle | Unbalanced/stray bold markers and broken sentence: "This means that any data you submit as a test will behave exactly **like a real submission, but since it's detached, nothing **will appear on the NCBI website". The ** markers are misplaced, creating a broken sentence. | Fix the bold markers, e.g. "...will behave exactly like a real submission, but since it's detached, **nothing will appear** on the NCBI website". |
| low | accessibility | Line 133, Customized Column Names tip | Link text "please refer to the column definitions in `task_submission.wdl`" points to a specific code line; acceptable, but "#L65" anchor and reading raw WDL is not friendly to the stated non-technical audience with no alternative provided. | Add a note that this is an advanced reference, or surface the key required/optional fields in the doc itself (partly done in the toggles below). |
| low | ease-of-use | Line 22, Terra_2_NCBI_PHB section | Acronyms NCBI, BioSample, SRA, FTP, XML are used; SRA and FTP in particular are unexpanded for non-technical readers (Sequence Read Archive; File Transfer Protocol). | Expand SRA and FTP on first use. |
| low | ease-of-use | Line 64, Collating BioSample Metadata section | MIxS acronym appears ("do not fit under the MIxS, Pathogen, or Virus packages") without expansion. | Expand MIxS (Minimum Information about any (x) Sequence) or link to a definition. |

### `docs/workflows/standalone/amr_search.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 14, AMR_Search_PHB section | AMR acronym is never expanded (antimicrobial resistance); the title and body use AMR throughout with no definition. | Expand on first use: "antimicrobial resistance (AMR) profiling". |
| low | grammar | Line 14, AMR_Search_PHB section | Missing article before AMRsearch: "utilizing `AMRsearch` tool from Pathogenwatch". | Insert "the": "utilizing the `AMRsearch` tool from Pathogenwatch". |
| low | grammar | Line 14, AMR_Search_PHB section | Typo in image alt text on line 11: "running paarsnp and parse_amr_json tools" - 'paarsnp' is likely a misspelling of PAARSNP/paarsnp; if it is the tool name it is fine, but the doubled 'a' looks like a typo. | Verify the tool name spelling; if it should be 'paarsnp' (single tool) confirm, otherwise correct the doubled 'a'. |
| low | accessibility | Lines 18-28, species/NCBI code table | The table is presented with minimal lead-in; line 16 mentions NCBI codes but the table's purpose (look up your species to find the code to enter) could be clearer for screen-reader users encountering the table. | Add one sentence directly before the table: "Find your species in the table below and use the corresponding NCBI Code as the workflow input." |

### `docs/workflows/standalone/cauris_cladetyper.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 9, Cauris_CladeTyper_PHB section | WGS acronym is used without expansion (whole genome sequencing); also "clade" is technical jargon for the non-technical audience. | Expand WGS on first use ("whole genome sequencing (WGS)") and add a brief plain-language note on what a clade is. |
| low | grammar | Line 9, Cauris_CladeTyper_PHB section | Redundant capitalization/phrasing: "The Cauris_CladeTyper_PHB Workflow is designed to assign the clade to..." - "assign the clade to" reads awkwardly. | Reword to "...is designed to assign a clade to Candidozyma auris ... WGS assemblies". |
| low | ease-of-use | Line 9, Cauris_CladeTyper_PHB section | The page lacks any "when/why to use" framing beyond one dense sentence, and there is no Inputs/Outputs context for a first-time user. | Add a short sentence on when a user would run this workflow (e.g. after assembling a C. auris genome) and what the result tells them. |

### `docs/workflows/standalone/concatenate_illumina_lanes.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 9, Concatenate_Illumina_Lanes_PHB section | "multi-lane FASTQ files" and "read type (forward or reverse)" assume familiarity with sequencing lanes and read directionality; brief context would help non-technical users know when this applies to their data. | Add a one-line note explaining that some sequencers split one sample's reads across multiple lane files, and this workflow merges them so downstream workflows can use them. |

### `docs/workflows/standalone/dorado_basecalling.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 9, Dorado_Basecalling_PHB section | "basecalling" and "POD5" are core concepts used immediately without definition; GPU is used without expansion. Non-technical readers may not know what basecalling is or why GPU acceleration matters. | Add a brief plain-language definition of basecalling (converting raw electrical signal into DNA sequence) and note POD5 is Oxford Nanopore's raw signal file format; expand GPU on first use. |
| low | grammar | Line 29, step 2 'Copy the GCS Path' | Sentence missing ending period: "...right-click the collection name and select \"Copy link address\"". | Add a closing period. |
| low | accessibility | Line 29, step 2 'Copy the GCS Path' | Instruction relies on a UI action with the screenshot, but the text "right-click the collection name and select 'Copy link address'" assumes the user knows what 'the collection name' refers to. | Clarify what 'the collection' is in the Terra Data Uploader (e.g. the named group/folder of uploaded files). |
| low | ease-of-use | Line 53, Model Type Selection section | "simplex model path or a model complex" and modification codes (5mCG_5hmCG, 6mA) are advanced jargon with no plain-language explanation of when a user would need them. | Add a note that most users should use sup/hac/fast and that manual model specification is an advanced option. |
| low | ease-of-use | Line 47, Model Type Selection | Inconsistent terminology: model types are referred to as both "model" and "model version" ("automatically select the appropriate model version"), which may confuse readers about whether they pick a model or a version. | Use consistent terminology, clarifying that the workflow auto-selects the best available version of the chosen model type. |

### `docs/workflows/standalone/gambit_query.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | accessibility | Line 31 | GAMBIT is expanded in the reference blockquote at the very bottom of the page, but the acronym is used much earlier (line 9, the page title, and section heading) without expansion at first use. | Expand GAMBIT (Genomic Approximation Method for Bacterial Identification and Tracking) at first use near line 9 so screen-reader users and newcomers get the definition in context. |
| low | ease-of-use | Line 9, GAMBIT_Query_PHB section | The single descriptive sentence assumes the reader knows what 'taxon assignment' and 'genome assembly' mean. There is no 'why/when to use this' guidance, unlike most other workflow pages which include a use-case note. | Add a short plain-language explanation, e.g. 'taxon assignment means identifying which species (or other taxonomic group) a sample belongs to' and a brief note on when a user would choose this standalone workflow. |

### `docs/workflows/standalone/kraken2.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | typo | Suggested databases table, 'viral w/human' row (line 41) | This row appears to be missing the Source column value: the cells are shifted so that '4.5' (size) lands in the Source column and the date/k-mer columns are off by one compared to every other row. The Source URL is missing. | Add the missing Source cell (or an empty cell) so the row's columns align with the header: Source, Database Size (GB), Date of Last Update, Bracken K-mer Lengths. |
| medium | ease-of-use | Suggested databases table, Date of Last Update column (lines 36-43) | Dates are written in ambiguous/inconsistent day-month-year format (e.g. '18/5/2022', '12/1/2024', '18/4/2023', '13/11/2020') while the rest of the docs (e.g. metabuli.md uses YYYY/MM/DD). A US public-health reader will likely misread '12/1/2024' as December 1 rather than 12 January. | Use an unambiguous, consistent date format such as ISO (2022-05-18) or spelled-out month (18 May 2022) across all rows. |
| low | typo | Line 93, Example Kraken2 report toggle | Spacing around the parenthetical em-dash phrase is inconsistent: 'very little -if any- read contamination' uses hyphens with surrounding spaces rather than proper dashes. | Use spaced en dashes or em dashes: 'very little — if any — read contamination'. |
| low | typo | Lines 93, 96, 105 etc. | Inconsistent spacing before percent sign: '84.35 %' and '6 %' (space) versus '~2 %' — and elsewhere in the repo style is typically no space ('98.78%'). | Standardize to no space before the percent sign (e.g. '84.35%') to match the rest of the documentation. |
| low | grammar | Suggested databases table, Kalamari row (line 36) | Comma splice / unnecessary comma: 'a database of complete public assemblies, that has been fine-tuned'. | Remove the comma before 'that': 'a database of complete public assemblies that has been fine-tuned for enteric pathogens'. |
| low | accessibility | Suggested databases table, Source column, Kalamari row (line 36) | The Source cell contains a bare '‣' glyph as the only content, which conveys no information to a screen reader and no meaning to a sighted reader either. | Replace the bare marker with a descriptive link or text (e.g. the Kalamari GitHub/source URL), or leave the cell empty if no source applies. |
| low | ease-of-use | Suggested databases table, two 'EuPathDB48' rows (lines 42-43) | Two different databases share the identical name 'EuPathDB48' (different dates/sizes: 2020-11-13 vs 2023-04-07). A non-technical reader cannot tell which to choose since the names are not distinguishable. | Differentiate the names, e.g. 'EuPathDB48 (2020)' and 'EuPathDB48 (2023)', so users can tell the two versions apart. |
| low | ease-of-use | Lines 11, 28 and elsewhere | Terms like 'k-mer', 'metagenomic', and 'horizontal gene transfer' are used without definition for a low-to-non-technical audience. | Add brief inline definitions on first use (e.g. 'k-mers, short DNA subsequences of length k') or link to a glossary. |

### `docs/workflows/standalone/metabuli.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | accessibility | Line 23 | Acronyms 'taxdump', 'RefSeq', and 'GTDB' are used without expansion on first use in the database-selection note. | Expand or briefly gloss these on first use (e.g. 'GTDB, the Genome Taxonomy Database') for a non-technical reader. |
| low | accessibility | Lines 109-111, Example Krona report toggle | The toggle says 'Below is an example of the krona_html for a bacterial sample', but the alt text and the example are the same generic Krona image reused from kraken2.md; the surrounding text describes a bacterial sample while the page's own example report (line 83) is an HIV viral sample, which may confuse readers about what they should be seeing. | Use a Metabuli-specific example or clarify in the text that this Krona image is a generic illustration of the report format, not the HIV example shown above. |
| low | ease-of-use | Line 11 | Dense technical sentence aimed at a non-technical audience: 'Metabuli uses a novel k-mer structure, called a "metamer", which incorporates both the DNA sequence for high specificity and amino acid conservation for sensitive homology detection.' Terms k-mer, specificity, amino acid conservation, and homology detection are undefined. | Either move this implementation detail to the Technical Details section, or add a plain-language summary first (e.g. 'Metabuli compares both DNA and protein information, which helps it identify organisms accurately even when sequences are not an exact match'). |
| low | ease-of-use | Line 13 | ONT is introduced as an acronym without expansion ('fastplong (ONT)'). A non-technical reader may not know ONT means Oxford Nanopore Technologies / long reads. | Expand on first use: 'fastplong (Oxford Nanopore Technologies, ONT, long reads)'. |

### `docs/workflows/standalone/ncbi_amrfinderplus.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 11 | Very long, dense paragraph packed with undefined jargon (Hidden Markov Models, reference gene catalog, docker image, workspace data element) for a non-technical reader. It mixes 'how it works' detail with practical advice in one block. | Break into shorter bullets/sentences and briefly define 'docker image' and 'workspace data element' (or link to the relevant Terra docs), since these are operational steps the user must take. |
| low | typo | Line 9 | Double space after the first sentence: 'stress genes.  Such AMR genes' and again 'point mutations.  In TheiaProk'. | Replace the double spaces with single spaces. |
| low | grammar | Line 11 | Informal phrasing 'You might like to save this docker image' is vague guidance for an instructional doc. | Rephrase as a clear recommendation, e.g. 'We recommend saving this docker image as a workspace data element to make future updates easier.' |
| low | accessibility | Line 11 | HMM is expanded once ('Hidden Markov Models (HMMs)') which is good, but 'AMR' is used in the title and line 9 before its expansion 'antimicrobial resistance (AMR)' which also appears on line 9 — verify AMR is expanded at true first use. | Ensure AMR is expanded the first time it appears in body text (line 9 'antimicrobial resistance (AMR) genes' is correct; confirm the heading/quick-facts usage does not precede a definition the reader can see). |

### `docs/workflows/standalone/ncbi_scrub.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 9 | Run-on / awkward sentence: 'NCBI Scrub ... is based on the SRA Taxonomy Analysis Tool that will take as input a FASTQ file, and produce as output a FASTQ file in which all reads identified as potentially of human origin are either removed (default) or masked with N.' The relative clause 'that will take as input...' grammatically attaches to the SRA Taxonomy Analysis Tool rather than to NCBI Scrub, and the comma before 'and produce' is unneeded. | Recast so the subject is NCBI Scrub, e.g. 'NCBI Scrub takes a FASTQ file as input and produces a FASTQ file in which all reads identified as potentially human are either removed (default) or masked with N.' Remove the comma before 'and'. |
| low | ease-of-use | Line 9-10 | No blank line between the end of the intro sentence (line 9) and 'There are two NCBI Scrub workflows:' (line 10), and the page gives no 'why/when to use this' context. Acronyms HRRT, SRA, and FASTQ are introduced without plain-language framing for the audience. | Add a blank line before 'There are two...' for clean rendering, and add a short sentence on when a user would run this (removing human reads before sharing or further analysis). Briefly gloss SRA and FASTQ. |
| low | ease-of-use | Line 32, Workflow Tasks | Uses the jargon verb 'dehost' ('one to dehost the input reads') which a non-technical reader may not recognize. | Use plain language or define it, e.g. 'one to remove human reads (dehost) the input reads and another to screen the clean reads with Kraken2 and the viral+human database.' |

### `docs/workflows/standalone/phylocompare.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 9, PhyloCompare_PHB section | Heavy phylogenetics jargon with no definitions for a non-technical audience: 'cophylogeny plot', 'tip arrangements', 'topologies', 'tip and branch arrangement', 'outgroup', 'midpoint'. The reader is not told why or when they would compare two trees. | Add a one-sentence plain-language purpose ('Use PhyloCompare to check whether two phylogenetic trees built from the same samples agree') and briefly define key terms (tip = a sample at the end of a branch; topology = the branching pattern; outgroup = a known distantly related sample used to root the tree). |
| low | grammar | Line 14 | Missing comma after the introductory clause: 'If no rooting options are supplied PhyloCompare will determine...' | Add a comma: 'If no rooting options are supplied, PhyloCompare will determine if the trees are rooted or unrooted.' |
| low | grammar | Line 11 | Comma splice joining two independent clauses: 'It is recommended to root a phylogeny and PhyloCompare can root upon an outgroup tip or the midpoint.' | Split into two sentences or add a comma before 'and': 'It is recommended to root a phylogeny, and PhyloCompare can root on an outgroup tip or the midpoint.' (also 'root upon' reads better as 'root on'). |
| low | ease-of-use | Line 19, phylovalidate_flag errors toggle | Very dense sentence relying on 'polytomy', 'non-0 length branches', and 'topology' without definition, likely impenetrable to the target audience. | Define 'polytomy' (a point where more than two branches split at once) and simplify the explanation of why it can confound distance calculations. |
| low | ease-of-use | Line 21 | Confusing wording: 'If flags are accompanied by a ">0" phylocompare_distance, then this indicates no distance was calculated.' A '>0' value normally means a distance WAS calculated, so this reads contradictorily. | Verify and clarify the intended condition (likely a sentinel/negative or 'NA' value indicates no distance), and reword so the relationship between the flag and the distance value is unambiguous. |

### `docs/workflows/standalone/rasusa.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 41, Verify Downsampling note | Awkward 'size/s' notation: 'comparing read file size/s to the original read file size/s'. | Write out plurals: 'comparing the read file size(s) to the original read file size(s)' or simply 'read file sizes'. |
| low | accessibility | Lines 15, 79 (RASUSA) — color/emoji-only reference | References to checking 'the box ✅ in the Terra workflow interface' rely on a green check emoji; the meaning depends on recognizing the icon. | Describe the control textually, e.g. 'the "Use call caching" checkbox in the Terra workflow interface', so the instruction is clear without relying on the emoji. |
| low | ease-of-use | Lines 9-12, RASUSA_PHB section | The workflow's core function (randomly subsampling/downsampling reads to a target coverage) is never stated up front; the page jumps straight to two use cases. The acronym LOD is expanded but 'coverage', 'subsample', and 'downsampling' are used without explanation. | Add a one-line description first, e.g. 'RASUSA randomly reduces (subsamples) the number of sequencing reads in a sample down to a target genome coverage.' Briefly define 'coverage'. |

### `docs/workflows/standalone/rename_fastq.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | ease-of-use | Line 9, Rename_FASTQ_PHB section | No 'why/when to use this' guidance. A non-technical reader is not told why they would need to rename FASTQ files (e.g. standardizing sample names before running other workflows). | Add a short sentence explaining a typical use case for renaming/recompressing FASTQ files. |

### `docs/workflows/standalone/tbprofiler_tngs.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | ease-of-use | Line 9, TBProfiler_tNGS_PHB section | The page gives no description of what the workflow actually does. It only states it is experimental. A reader cannot tell what TBProfiler_tNGS is for, what tNGS means, or what inputs/outputs to expect. tNGS (targeted Next-Generation Sequencing) and TBProfiler are unexpanded. | Even for an experimental workflow, add one or two sentences describing its purpose (e.g. profiling Mycobacterium tuberculosis drug resistance from targeted NGS data) and expand tNGS on first use. |

### `docs/workflows/standalone/theiavalidate.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Line 38 | Subject-verb agreement: 'values will pass if at least one criteria is met' — 'criteria' is plural; the singular is 'criterion'. | Change to 'if at least one criterion is met'. |
| low | accessibility | Lines 25-29, Validation Starting Points table | The 'Columns to Compare' cells contain extremely long, unbroken comma-delimited strings of dozens of column names with no spaces. This is very hard to read for everyone and especially difficult for screen-reader users, and the table provides no textual explanation of what these column lists represent. | Move the long column lists to the linked TSV files (already provided) or format them as wrapped/bulleted lists, and add a sentence explaining these are example metric sets Theiagen uses. |
| low | accessibility | Lines 126-134, Example Data and Outputs | The output example links use the bare file name as link text (e.g. 'example_summary.pdf'), which is acceptable, but the introductory line 'please observe the following example and outputs:' and the list give no description of what each file contains, so a screen-reader user gets a list of filenames with no context. | Add a short description after each link (e.g. 'example_summary.pdf — the human-readable PDF report summarizing all differences'). |
| low | ease-of-use | Line 15, TheiaValidate_PHB section | States 'A summary PDF report is produced in addition to an Excel spreadsheet', but the rest of the page and the diagram (line 13) and outputs refer to TSV files, not Excel. This is inconsistent terminology for the same output and may confuse users about the file format. | Use consistent terminology — clarify whether the tabular output is a TSV or an Excel/.xlsx file, and use the same term throughout. |
| low | ease-of-use | Line 36, PERCENT_DIFF bullet | Potentially confusing example: 'if the decimal percentage is 0.02, the test will indicate a failure if the values ... are more than 2% different.' The relationship between the decimal 0.02 and '2%' may not be obvious to a non-technical reader. | Clarify that 0.02 is expressed as a fraction equal to 2% (e.g. '0.02 = 2%'), so the mapping is explicit. |

### `docs/workflows_overview/workflows_alphabetically.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| low | grammar | Lines 16-17 | Inconsistent terminal punctuation in the abbreviation/definition list: 'Sample-level' and 'Set-level' definitions end without a period, while 'Table-level' ends with a period. | Make terminal punctuation consistent across the three definitions (all with or all without a period). |

### `docs/workflows_overview/workflows_kingdom.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | accessibility | Lines 9, 17, 25, 33 (section headings) | Heading hierarchy skips a level: the page title is H1 (frontmatter) and the first content heading jumps to H3 (### Any taxa) with no H2 in between. Skipped heading levels hinder screen-reader navigation. | Use H2 (##) for the top-level category sections, or add an intervening H2, so the hierarchy does not jump from H1 to H3. |
| low | grammar | Lines 42-43 | Inconsistent terminal punctuation in the level definitions ('Sample-level'/'Set-level' have no period; 'Table-level' has one). | Standardize punctuation across the definitions. |

### `docs/workflows_overview/workflows_type.md`

| Sev | Category | Location | Issue | Suggestion |
| --- | --- | --- | --- | --- |
| medium | accessibility | Lines 9, 17, 25, 33, 41, 49, 57, 65 (section headings) | Heading hierarchy skips a level: page title is H1 (frontmatter) and content headings jump directly to H3 (### Data Import, etc.) with no H2. This breaks logical heading nesting for assistive technology. | Promote the category sections to H2, or insert an H2 above them, to avoid the H1-to-H3 jump. |
| low | grammar | Lines 74-75 | Inconsistent terminal punctuation in the workflow-level definitions ('Sample-level'/'Set-level' lack a period; 'Table-level' has one). | Standardize punctuation across the definitions. |
