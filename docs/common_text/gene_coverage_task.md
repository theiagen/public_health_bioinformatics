---
title: Task Fragment `gene_coverage`
fragment: true
---
??? task "`gene_coverage`: Depth and Breadth of Coverage Calculations"
    This task calculates average read depth and the percent of a region (typically genes) covered above a minimum depth and quality using samtools (Pysam) and basic arithmetic.

    Outputs are reported as JSON-based "maps" that relate gene names to their breadth and depth of coverage. For example:

    `percent_coverage_by_gene`:
    ```
    {
        "GENE1": 99.5
        "GENE2": 10
    }
    ```

    `depth_by_gene`:
    ```
    {
        "GENE1": 10,
        "GENE2": 1
    }

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
    ??? dna "GBFF and BED file usage"
        For fungal pathogens, either a GBFF or a BED file may be used for gene coverage coordinate selection.

        - If a GBFF is used, a comma-delimited `query_genes` list must be provided (for example: `geneA,geneB,geneC`) to extract gene coordinates
        - If a BED is used, gene names can be taken from the BED entries; if `query_genes` is supplied, particular regions will be extracted

        The following query genes are used by default:

        - *Aspergillus fumigatus*: `Cyp51A`, `HapE`, `AFUA_4G08340` (COX10 in the default reference)
        - *Candidozyma auris*: `FKS1`, `lanosterol.14-alpha.demethylase`, `uracil.phosphoribosyltransferase`, `B9J08_005340`, `B9J08_000401`, `B9J08_003102`, `B9J08_003737`, `B9J08_005343`
        - *Cryptococcus neoformans*: `CNA00300` (ERG11 in the default reference)
<!-- endif -->

    !!! techdetails "Gene Coverage Technical Details"        
        |  | Links |
        | --- | --- |
        | Task | [task_gene_coverage.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_gene_coverage.wdl) |
        | Software Source Code | [samtools on GitHub](https://github.com/samtools/samtools) |
        | Software Documentation | [samtools Manual](https://www.htslib.org/doc/samtools.html) |
        | Original Publication(s) | [Twelve years of SAMtools and BCFtools](https://doi.org/10.1093/gigascience/giab008) |
