??? task "Snippy_Variants"
    ##### Snippy_Variants
<!-- if: snippy_streamline -->
    `Snippy_Variants` uses Snippy to align the assemblies for each sample against the reference genome to call SNPs, MNPs and INDELs according to optional input parameters. 
<!-- endif -->

<!-- if: snippy_variants -->
    `Snippy_Variants` uses Snippy to align reads to the reference and call SNPs, MNPs and INDELs according to optional input parameters.
<!-- endif -->

    Optionally, if the user provides a value for `query_gene`, the variant file will be searched for any mutations in the specified regions or annotations. The query string MUST match the gene name or annotation as specified in the GenBank file and provided in the output variant file in the `snippy_results` column.

    ??? toggle "QC Metrics from Snippy_Variants"
<!-- if: snippy_streamline -->
        !!! warning 
            The following QC metrics may not be applicable to your dataset as they are geared towards read data, not assemblies. Use these metrics with caution.
<!-- endif -->

        This task also extracts QC metrics from the Snippy output for each sample and saves them in per-sample TSV files (`snippy_variants_qc_metrics`). These per-sample QC metrics include the following columns:

        - **samplename**: The name of the sample.
        - **reads_aligned_to_reference**: The number of reads that aligned to the reference genome.
        - **total_reads**: The total number of reads in the sample.
        - **percent_reads_aligned**: The percentage of reads that aligned to the reference genome.
        - **variants_total**: The total number of variants detected between the sample and the reference genome.
        - **percent_ref_coverage**: The percentage of the reference genome covered by reads with a depth greater than or equal to the `min_coverage` threshold (default is 10).
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

<!-- if: snippy_variants -->
    !!! tip "QC Metrics for Phylogenetic Analysis"
        These QC metrics provide valuable insights into the quality and coverage of your sequencing data relative to the reference genome. Monitoring these metrics can help identify samples with low coverage, poor alignment, or potential issues that may affect downstream analyses, and we recommend examining them before proceeding with phylogenetic analysis if performing Snippy_Variants and Snippy_Tree separately.

        These per-sample QC metrics can also be combined into a single file (`snippy_combined_qc_metrics`) in downstream workflows, such as `snippy_tree`, providing an overview of QC metrics across all samples.
<!-- endif -->

    !!! techdetails "Snippy Variants Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_snippy_variants.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snippy_variants.wdl)<br>[task_snippy_gene_query.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snippy_gene_query.wdl) |
        | Software Source Code | [Snippy on GitHub](https://github.com/tseemann/snippy) |
        | Software Documentation | [Snippy on GitHub](https://github.com/tseemann/snippy) |
