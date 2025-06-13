??? task "`metabuli`"

    The `metabuli` task is used to classify and extract reads against a reference database. Metabuli uses a novel k-mer structure, called metamer, to analyze both amino acid (AA) and DNA sequences. It leverages AA conservation for sensitive homology detection and DNA mutations for specific differentiation between closely related taxa.

    ??? `cpus` / `memory`
        Increasing the memory and cores allocated to Metabuli can substantially increase throughput.

<!-- if: theiaviral -->
    ??? dna "`extract_unclassified`"
        This parameter determines whether unclassified reads should also be extracted and combined with the `taxon`-specific extracted reads. By default, this is set to `false`, meaning that only reads classified to the specified input `taxon` will be extracted.

    ???+ warning "Important"
        This task will extract reads classified to the input `taxon` and **all of its descendant taxa**. The `rank` input parameter controls the extraction of reads classified at the specified `rank` and all suboridante taxonomic levels. See task `ncbi_identify` under the **Taxonomic Identification** section for more details on the `rank` input parameter.
<!-- endif -->

    !!! techdetails "Metabuli Technical Details"
    |  | Links |
    | --- | --- |
    | Task | [task_metabuli.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/contamination/task_metabuli.wdl) |
    | Software Source Code | [Metabuli on GitHub](https://github.com/steineggerlab/Metabuli) |
    | Software Documentation | [Metabuli Documentation](https://github.com/steineggerlab/Metabuli/blob/master/README.md) |
    | Original Publication(s) | [Metabuli Paper](https://doi.org/10.1038/s41592-024-02273-y) |