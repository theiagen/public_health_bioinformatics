??? task "`metabuli`"

    The `metabuli` task is used to classify and optionally extract reads against a reference database. Metabuli uses a novel k-mer structure, called metamer, to analyze both amino acid (AA) and DNA sequences. It leverages AA conservation for sensitive homology detection and DNA mutations for specific differentiation between closely related taxa.

<!-- if: metabuli -->
    ??? dna "`taxon_id` input parameter"
        `taxon_id` triggers read extraction by retrieving the inputted NCBI taxon ID and all descendant taxon IDs derived from the input.

    ??? dna "Precision mode and `min_score` / `min_sp_score` input parameters"
        The `min_score` parameter is the minimum score (DNA-level identity) required for a read to be classified and the `min_sp_score` parameter is the minimum score for a read to be classified at or below species rank. Metabuli precision mode is defined by its authors as more stringently setting the `min_score` and `min_sp_score` parameters for specific read types:

        - Illumina short reads: `min_score` = 0.15, `min_sp_score` = 0.5
        - ONT long reads: `min_score` = 0.008

<!-- endif -->

    ??? dna "`taxdump_path` input parameter"
        The `taxdump_path` directs the task toward a taxonkit-generated taxdump file, e.g. [from NCBI](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/) or [from GTDB](https://github.com/shenwei356/gtdb-taxdump/releases). This is not necessary to edit unless users want a more recent taxdump than what Theiagen hosts, or if users want to reference a different taxonomy. By default, Theiagen uses the NCBI taxonomy hierarchy.

    ??? dna "`cpu` / `memory` input parameters"
        Increasing the memory and cpus allocated to Metabuli can substantially increase throughput.

    ??? dna "`extract_unclassified` input parameter"
        This parameter determines whether unclassified reads should also be extracted and combined with the `taxon`-specific extracted reads. By default, this is set to `false`, meaning that only reads classified to the specified input `taxon` will be extracted.

<!--if: theiaviral -->
    ???+ warning "Descendant taxa reads are extracted"
        This task will extract reads classified to the input `taxon` and **all of its descendant taxa**. The `rank` input parameter controls the extraction of reads classified at the specified `rank` and all subordiante taxonomic levels. See task `ncbi_identify` under the **Taxonomic Identification** section above for more details on the `rank` input parameter.
<!-- endif -->

    !!! techdetails "Metabuli Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_metabuli.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/contamination/task_metabuli.wdl) |
        | Software Source Code | [Metabuli on GitHub](https://github.com/steineggerlab/Metabuli) |
        | Software Documentation | [Metabuli Documentation](https://github.com/steineggerlab/Metabuli/blob/master/README.md) |
        | Original Publication(s) | [Metabuli: sensitive and specific metagenomic classification via joint analysis of amino acid and DNA](https://doi.org/10.1038/s41592-024-02273-y) |