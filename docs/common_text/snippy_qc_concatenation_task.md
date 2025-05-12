
??? task "Snippy_Variants QC Metrics Concatenation (optional)"
    ##### Snippy_Variants QC Metric Concatenation (optional)

    Optionally, the user can provide the `snippy_variants_qc_metrics` file produced by the Snippy_Variants workflow as input to the workflow to concatenate the reports for each sample in the tree. These per-sample QC metrics include the following columns:

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

    The combined QC metrics file includes the same columns as above for all samples. Note that the last set of columns (`#rname` to `meanmapq`) may repeat for each chromosome or contig in the reference genome.

    !!! tip "QC Metrics for Phylogenetic Analysis"
        These QC metrics provide valuable insights into the quality and coverage of your sequencing data relative to the reference genome. Monitoring these metrics can help identify samples with low coverage, poor alignment, or potential issues that may affect downstream analyses, and we recommend examining them before proceeding with phylogenetic analysis if performing Snippy_Variants and Snippy_Tree separately.

    !!! techdetails "Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_cat_files.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/file_handling/task_cat_files.wdl) |
