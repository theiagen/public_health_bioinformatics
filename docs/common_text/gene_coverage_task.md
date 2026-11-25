---
title: Task Fragment `gene_coverage`
fragment: true
---
??? task "`gene_coverage`: Depth and Breadth of Coverage Calculations"
    This task calculates average read depth, the percent of a region (typically genes) covered above a minimum depth and quality (breadth), and number of reads mapped using samtools (Pysam) and basic arithmetic.

    Outputs are reported as JSON-based "maps" that relate gene names to the statistic of interest. For example:

    `depth_by_gene`:
    ```
    {
        "GENE1": 10,
        "GENE2": 1
    }
    ```

<!-- if: viral -->
    By default, this task runs for SARS-CoV-2 and Mpox.

    !!! warning "Region coordinates must be relevant to the reference genome"
        Please note that default BEDfiles contain gene coordinates that may not directly match user-provided or dynamically-selected reference genomes (TheiaViral).

    ??? dna "BED file usage"
        In viral characterization workflows, gene coverage regions are supplied with a BED file.

        - To extract custom regions of interest, populate the `reference_gene_locations_bed` input (task `theiacov` / `morgana_magic`)
        - If no custom BED is provided, organism defaults are used when available
        - BED files should include a gene name in column 4 to label output
<!-- endif -->

<!-- if: theiaeuk -->
    ??? dna "GFF and BED file usage"
        Gene coverage runs when a reference GFF (default or user-supplied) is available alongside either `query_genes` or `query_genes_bed`, which are mutually exclusive:

        - **`query_genes`**: the comma-delimited gene names (for example: `geneA,geneB,geneC`) are matched against the product and locus tag fields of the GFF, and the CDS coordinates of each matched gene are quantified
        - **`query_genes_bed`**: takes precedence over `query_genes` and any organism default; each BED region is quantified as written and labelled by the BED name column (column 4), so the GFF does not contribute coordinates
<!-- endif -->

    ??? dna "How is coverage quantified?"
        Depth is quantified as the per-base average number of reads that primarily align to a given region and meet minimum quality thresholds. Breadth of coverage (percent reference coverage) is quantified as the percent of a region with bases that are covered with a depth that meets or exceeds the minimum (10 by default).

    !!! techdetails "Gene Coverage Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_gene_coverage.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_gene_coverage.wdl) |
        | Software Source Code | [theiagene on GitHub](https://github.com/theiagen/theiagene) |
        | Software Documentation | [samtools Manual](https://www.htslib.org/doc/samtools.html) |
        | Original Publication(s) | [Twelve years of SAMtools and BCFtools](https://doi.org/10.1093/gigascience/giab008) |
