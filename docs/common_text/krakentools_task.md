??? task "`krakentools`"

    The `task_krakentools.wdl` task extracts reads from the Kraken2 output file. It uses the [KrakenTools](https://github.com/jenniferlu717/KrakenTools) package to extract reads classified at any user-specified taxon ID.

<!-- if: theiaviral -->
    ??? dna "`extract_unclassified`"
        This parameter determines whether unclassified reads should also be extracted and combined with the `taxon`-specific extracted reads. By default, this is set to `false`, meaning that only reads classified to the specified input `taxon` will be extracted.

    ???+ warning "Important"
        This task will extract reads classified to the input `taxon` and **all of its descendant taxa**. The `rank` input parameter controls the extraction of reads classified at the specified `rank` and all suboridante taxonomic levels. See task `ncbi_identify` under the **Taxonomic Identification** section for more details on the `rank` input parameter.
<!-- endif -->

    !!! techdetails "KrakenTools Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_krakentools.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/task_krakentools.wdl) |
        | Software Source Code | [KrakenTools on GitHub](https://github.com/jenniferlu717/KrakenTools) |
        | Software Documentation | [KrakenTools](https://github.com/jenniferlu717/KrakenTools/blob/master/README.md) |
        | Original Publication(s) | [Metagenome analysis using the Kraken software suite](https://doi.org/10.1126/scitranslmed.aap9489) |