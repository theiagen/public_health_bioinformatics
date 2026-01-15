<!-- if: snippy_streamline|snippy_variants -->
??? task "Snippy_Variants"
    ##### Snippy_Variants
<!-- endif -->
<!-- if: snippy_streamline -->
    `Snippy_Variants` uses Snippy to align the assemblies for each sample against the reference genome to call SNPs, MNPs and INDELs according to optional input parameters. 

<!-- endif -->
<!-- if: snippy_variants -->
    `Snippy_Variants` uses Snippy to align reads to the reference and call SNPs, MNPs and INDELs according to optional input parameters.

<!-- endif -->
<!-- if: snippy_streamline|snippy_variants -->
    Optionally, if the user provides a value for `query_gene`, the variant file will be searched for any mutations in the specified regions or annotations. The query string MUST match the gene name or annotation as specified in the GenBank file and provided in the output variant file in the `snippy_results` column.

    ??? toggle "QC Metrics from Snippy_Variants"
<!-- endif -->
<!-- if: snippy_streamline -->
        !!! warning 
            The following QC metrics may not be applicable to your dataset as they are geared towards read data, not assemblies. Use these metrics with caution.

<!-- endif -->
<!-- if: snippy_streamline|snippy_variants -->
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

<!-- endif -->
<!-- if: snippy_variants -->
    !!! tip "QC Metrics for Phylogenetic Analysis"
        These QC metrics provide valuable insights into the quality and coverage of your sequencing data relative to the reference genome. Monitoring these metrics can help identify samples with low coverage, poor alignment, or potential issues that may affect downstream analyses, and we recommend examining them before proceeding with phylogenetic analysis if performing Snippy_Variants and Snippy_Tree separately.

        These per-sample QC metrics can also be combined into a single file (`snippy_combined_qc_metrics`) in downstream workflows, such as `snippy_tree`, providing an overview of QC metrics across all samples.

<!-- endif -->
<!-- if: cauris|calbicans|afumigatus|cneoformans -->
??? task "`Snippy_Variants`: Antifungal Resistance Detection"
    To detect mutations that may confer antifungal resistance, `Snippy` is used to find all variants relative to the clade-specific reference, then these variants are queried for product names associated with resistance. It's important to note that unlike `amr_search`, this task reports all variants found in the searched targets.
<!-- endif -->

<!-- if: cauris -->
    - FKS1
    - ERG11 (lanosterol 14-alpha demethylase)
    - FUR1 (uracil phosphoribosyltransferase)
<!-- endif -->
<!-- if: calbicans -->
    - ERG11
    - GCS1 (FKS1)
    - FUR1
    - RTA2
<!-- endif -->
<!-- if: afumigatus -->
    - Cyp51A
    - HapE
    - COX10 (AFUA_4G08340)
<!-- endif -->
<!-- if: cneoformans -->
    - ERG11 (CNA00300)
<!-- endif -->

<!-- if: cauris|calbicans|afumigatus|cneoformans -->
    We query `Snippy` results to see if any mutations were identified in those genes. By default, we automatically check for the following loci (which can be overwritten by the user). You will find the mutations next to the locus tag in the `theiaeuk_snippy_variants_hits` column corresponding gene name (see below):
<!-- endif -->

<!-- if: cauris -->
    | **TheiaEuk Search Term** | **Corresponding Gene Name** |
    |---|---|
    | B9J08_005340 | ERG6 |
    | B9J08_000401 | FLO8 |
    | B9J08_005343 | Hypothetical protein (PSK74852) |
    | B9J08_003102 | MEC3 |
    | B9J08_003737 | ERG3 |
    | lanosterol.14-alpha.demethylase | ERG11 |
    | uracil.phosphoribosyltransferase | FUR1 |
    | FKS1 | FKS1 |    

    ??? toggle "Known resistance-conferring mutations for _Candidozyma auris_"
        Mutations in these genes that are known to confer resistance are shown below

        | **Organism** | **Found in** | **Gene name** | **Gene locus** | **AA mutation** | **Drug** | **Reference** |
        | --- | --- | --- | --- | --- | --- | --- |
        | Candidozyma auris | Human | ERG11 |  | Y132F | Fluconazole | [Simultaneous Emergence of Multidrug-Resistant _Candida auris_ on 3 Continents Confirmed by Whole-Genome Sequencing and Epidemiological Analyses](https://academic.oup.com/cid/article/64/2/134/2706620/Simultaneous-Emergence-of-Multidrug-Resistant) |
        | Candidozyma auris | Human | ERG11 |  | K143R | Fluconazole | [Simultaneous Emergence of Multidrug-Resistant _Candida auris_ on 3 Continents Confirmed by Whole-Genome Sequencing and Epidemiological Analyses](https://academic.oup.com/cid/article/64/2/134/2706620/Simultaneous-Emergence-of-Multidrug-Resistant) |
        | Candidozyma auris | Human | ERG11 |  | F126T | Fluconazole | [Simultaneous Emergence of Multidrug-Resistant _Candida auris_ on 3 Continents Confirmed by Whole-Genome Sequencing and Epidemiological Analyses](https://academic.oup.com/cid/article/64/2/134/2706620/Simultaneous-Emergence-of-Multidrug-Resistant) |
        | Candidozyma auris | Human | FKS1 |  | S639P | Micafungin  | [Activity of CD101, a long-acting echinocandin, against clinical isolates of Candida auris](https://www.sciencedirect.com/science/article/pii/S0732889317303498) |
        | Candidozyma auris | Human | FKS1 |  | S639P | Caspofungin | [Activity of CD101, a long-acting echinocandin, against clinical isolates of Candida auris](https://www.sciencedirect.com/science/article/pii/S0732889317303498) |
        | Candidozyma auris | Human | FKS1 |  | S639P | Anidulafungin | [Activity of CD101, a long-acting echinocandin, against clinical isolates of Candida auris](https://www.sciencedirect.com/science/article/pii/S0732889317303498) |
        | Candidozyma auris | Human | FKS1 |  | S639F | Micafungin | [A multicentre study of antifungal susceptibility patterns among 350 _Candida auris_ isolates (2009–17) in India: role of the ERG11 and FKS1 genes in azole and echinocandin resistance](https://academic.oup.com/jac/advance-article/doi/10.1093/jac/dkx480/4794718) |
        | Candidozyma auris | Human | FKS1 |  | S639F | Caspofungin | [A multicentre study of antifungal susceptibility patterns among 350 _Candida auris_ isolates (2009–17) in India: role of the ERG11 and FKS1 genes in azole and echinocandin resistance](https://academic.oup.com/jac/advance-article/doi/10.1093/jac/dkx480/4794718) |
        | Candidozyma auris | Human | FKS1 |  | S639F | Anidulafungin | [A multicentre study of antifungal susceptibility patterns among 350 _Candida auris_ isolates (2009–17) in India: role of the ERG11 and FKS1 genes in azole and echinocandin resistance](https://academic.oup.com/jac/advance-article/doi/10.1093/jac/dkx480/4794718) |
        | Candidozyma auris | Human | FUR1 | CAMJ_004922 | F211I | 5-flucytosine | [Genomic epidemiology of the UK outbreak of the emerging human fungal pathogen Candida auris](https://doi.org/10.1038/s41426-018-0045-x) |
<!-- endif -->
<!-- if: calbicans -->
    | **TheiaEuk Search Term** | **Corresponding Gene Name** |
    |---|---|
    | ERG11 | ERG11 |
    | GCS1 | FKS1 |
    | FUR1 | FUR1 |
    | RTA2 | RTA2 |
<!-- endif -->
<!-- if: afumigatus -->
    | **TheiaEuk Search Term** | **Corresponding Gene Name** |
    |---|---|
    | Cyp51A | Cyp51A |
    | HapE | HapE |
    | AFUA_4G08340 | COX10 |
<!-- endif -->
<!-- if: cneoformans -->
    | **TheiaEuk Search Term** | **Corresponding Gene Name** |
    |---|---|
    | CNA00300 | ERG11 |
<!-- endif -->

<!-- if: cauris|calbicans|afumigatus|cneoformans -->
    ??? toggle "Example Output Interpretation"
        For example, one sample may have the following output for the `theiaeuk_snippy_variants_hits` column:

        ```plaintext
        lanosterol.14-alpha.demethylase: lanosterol 14-alpha demethylase (missense_variant c.428A>G p.Lys143Arg; C:266 T:0),B9J08_000401: hypothetical protein (stop_gained c.424C>T p.Gln142*; A:70 G:0)
        ```

        Based on this, we can tell that ERG11 has a missense variant at position 143 (Lysine to Arginine) and B9J08_000401 (which is FLO8) has a stop-gained variant at position 142 (Glutamine to Stop).
<!-- endif -->

    !!! techdetails "Snippy Variants Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_snippy_variants.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snippy_variants.wdl)<br>[task_snippy_gene_query.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snippy_gene_query.wdl) |
        | Software Source Code | [Snippy on GitHub](https://github.com/tseemann/snippy) |
        | Software Documentation | [Snippy on GitHub](https://github.com/tseemann/snippy) |
