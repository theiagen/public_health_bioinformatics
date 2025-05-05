# Snippy_Variants

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Phylogenetic Construction](../../workflows_overview/workflows_type.md/#phylogenetic-construction) | [Bacteria](../../workflows_overview/workflows_kingdom.md/#bacteria), [Mycotics](../../workflows_overview/workflows_kingdom.md#mycotics), [Viral](../../workflows_overview/workflows_kingdom.md/#viral) | PHB v2.3.0 | Yes | Sample-level |

## Snippy_Variants_PHB

The `Snippy_Variants` workflow aligns single-end or paired-end reads (in FASTQ format), or assembled sequences (in FASTA format), against a reference genome, then identifies single-nucleotide polymorphisms (SNPs), multi-nucleotide polymorphisms (MNPs), and insertions/deletions (INDELs) across the alignment. If a GenBank file is used as the reference, mutations associated with user-specified query strings (e.g. genes of interest) can additionally be reported to the Terra data table.

!!! caption "Snippy_Variants Workflow Diagram"
    ![Snippy_Variants Workflow Diagram](../../assets/figures/Snippy_Variants.png)

!!! tip "Example Use Cases"
    - **Finding mutations** (SNPs, MNPs, and INDELs) in your own sample's reads relative to a reference, e.g. mutations in genes of phenotypic interest.
    - **Quality control:** When undertaking quality control of sequenced isolates, it is difficult to identify contamination between multiple closely related genomes using the conventional approaches in TheiaProk (e.g. isolates from an outbreak or transmission cluster). Such contamination may be identified as allele heterogeneity at a significant number of genome positions. `Snippy_Variants` may be used to identify these heterogeneous positions by aligning reads to the assembly of the same reads, or to a closely related reference genome and lowering the thresholds to call SNPs.
    - **Assessing support for a mutation**: `Snippy_Variants` produces a BAM file of the reads aligned to the reference genome. This BAM file can be visualized in IGV (see Theiagen Office Hours recordings) to assess the position of a mutation in supporting reads, or if the assembly of the reads was used as a reference, the position in the contig.
        - Mutations that are only found at the ends of supporting reads may be an error of sequencing.
        - Mutations found at the end of contigs may be assembly errors.

### Inputs

- Single or paired-end reads resulting from Illumina or IonTorrent sequencing can be used. For single-end data, simply omit a value for `read2`
- Assembled genomes can be used. Use the `assembly_fasta` input and omit `read1` and `read2`
- The reference file should be in fasta (e.g. `.fa`, `.fasta`) or [full GenBank](https://github.com/tseemann/snippy/issues/463#issuecomment-863344618) (`.gbk`) format. The mutations identified by Snippy_Variants are highly dependent on the choice of reference genome. Mutations cannot be identified in genomic regions that are present in your query sequence and not the reference.

!!! info "Query String"
    The query string can be a gene or any other annotation that matches the GenBank file/output VCF **EXACTLY**

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/input_tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="Snippy_Variants", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Workflow Tasks

`Snippy_Variants` uses Snippy to align reads to the reference and call SNPs, MNPs and INDELs according to optional input parameters. The output includes a file of variants that is then queried using the `grep` bash command to identify any mutations in specified genes or annotations of interest. The query string MUST match the gene name or annotation as specified in the GenBank file and provided in the output variant file in the `snippy_results` column.

!!! info "Quality Control Metrics"
    Additionally, `Snippy_Variants` extracts quality control (QC) metrics from the Snippy output for each sample. These per-sample QC metrics are saved in TSV files (`snippy_variants_qc_metrics`). The QC metrics include:

    - **samplename**: The name of the sample.
    - **reads_aligned_to_reference**: The number of reads that aligned to the reference genome.
    - **total_reads**: The total number of reads in the sample.
    - **percent_reads_aligned**: The percentage of reads that aligned to the reference genome; also available in the `snippy_variants_percent_reads_aligned` output column.
    - **variants_total**: The total number of variants detected between the sample and the reference genome.
    - **percent_ref_coverage**: The percentage of the reference genome covered by reads with a depth greater than or equal to the `min_coverage` threshold (default is 10); also available in the `snippy_variants_percent_ref_coverage` output column.
    - **#rname**: Reference sequence name (e.g., chromosome or contig name).
    - **startpos**: Starting position of the reference sequence.
    - **endpos**: Ending position of the reference sequence.
    - **numreads**: Number of reads covering the reference sequence.
    - **covbases**: Number of bases with coverage.
    - **coverage**: Percentage of the reference sequence covered (depth ≥ 1).
    - **meandepth**: Mean depth of coverage over the reference sequence.
    - **meanbaseq**: Mean base quality over the reference sequence.
    - **meanmapq**: Mean mapping quality over the reference sequence.
  
    Note that the last set of columns (`#rname` to `meanmapq`) may repeat for each chromosome or contig in the reference genome.

!!! tip "QC Metrics for Phylogenetic Analysis"
    These QC metrics provide valuable insights into the quality and coverage of your sequencing data relative to the reference genome. Monitoring these metrics can help identify samples with low coverage, poor alignment, or potential issues that may affect downstream analyses, and we recommend examining them before proceeding with phylogenetic analysis if performing Snippy_Variants and Snippy_Tree separately.

    These per-sample QC metrics can also be combined into a single file (`snippy_combined_qc_metrics`) in downstream workflows, such as `snippy_tree`, providing an overview of QC metrics across all samples.

!!! techdetails "Snippy Variants Technical Details"
    |  | Links |
    | --- | --- |
    | Task | [task_snippy_variants.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snippy_variants.wdl)<br>[task_snippy_gene_query.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snippy_gene_query.wdl) |
    | Software Source Code | [Snippy on GitHub](https://github.com/tseemann/snippy) |
    | Software Documentation | [Snippy on GitHub](https://github.com/tseemann/snippy) |

### Outputs

!!! tip "Visualize your outputs in IGV"
    Output bam/bai files may be visualized using IGV to manually assess read placement and SNP support.

!!! warning "Note on coverage calculations"
    The outputs from `samtools coverage` (found in the `snippy_variants_coverage_tsv` file) may differ from the `snippy_variants_percent_ref_coverage` due to different calculation methods. `samtools coverage` computes genome-wide coverage metrics (e.g., the proportion of bases covered at depth ≥ 1), while `snippy_variants_percent_ref_coverage` uses a user-defined minimum coverage threshold (default is 10), calculating the proportion of the reference genome with a depth greater than or equal to this threshold.

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/tables/all_outputs.tsv", input_table=False, filter_column="Workflow", filter_values="Snippy_Variants", columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///

