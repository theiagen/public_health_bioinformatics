<!-- if: theiaprok -->
??? task "`Kraken2`: Read Identification (optional)"
<!-- endif -->
<!-- if: theiacov|freyja -->
??? task "`Kraken2`: Read Identification"
<!-- endif -->
<!-- if: kraken -->
??? task "Kraken2"
<!-- endif -->
<!-- if: theiaviral -->
??? toggle "Read Identification"
<!-- endif -->

    `Kraken2` is a bioinformatics tool originally designed for metagenomic applications. It has additionally proven valuable for validating taxonomic assignments and checking contamination of single-species (e.g. bacterial isolate, eukaryotic isolate, viral isolate, etc.) whole genome sequence data.

<!-- if: theiacov|freyja -->
    Kraken2 is run on both the raw and clean reads.
<!-- endif -->

<!-- if: theiaviral -->
    This task runs on cleaned reads passed from the `read_QC_trim` subworkflow and outputs a Kraken2 report detailing taxonomic classifications. It also separates classified reads from unclassified ones.
<!-- endif -->

    !!! info "Database-dependent"
        This workflow automatically uses a viral-specific Kraken2 database. This database was generated in-house from RefSeq's viral sequence collection and human genome GRCh38. It's available at `gs://theiagen-public-resources-rp/reference_data/databases/kraken2/kraken2_humanGRCh38_viralRefSeq_20240828.tar.gz`.
<!-- endif -->
  
<!-- if: theiaprok -->
    As an alternative to `MIDAS` (see above), the `Kraken2` task can also be turned on through setting the `call_kraken` input variable as `true` for the identification of reads to detect contamination with non-target taxa.

    A database must be provided if this optional module is activated, through the kraken_db optional input. A list of suggested databases can be found on [Kraken2 standalone documentation](../standalone/kraken2.md#databases).
<!-- endif -->

<!-- if: kraken -->
    This workflow is database dependent, and one is required to run this task. Please see above for a list of suggested databases to provide through the `kraken2_db` input variable.
<!-- endif -->

<!-- if: theiameta -->
    Kraken2 is run on the set of raw reads, provided as input, as well as the set of clean reads that are resulted from the `read_QC_trim` workflow

    The Kraken2 software is database-dependent and **taxonomic assignments are highly sensitive to the database used**. An appropriate database should contain the expected organism(s) (e.g. _Escherichia coli_) and other taxa that may be present in the reads (e.g. _Citrobacter freundii_, a common contaminant).
<!-- endif -->

    !!! techdetails "Kraken2 Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_kraken2.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/contamination/task_kraken2.wdl) |
        | Software Source Code | [Kraken2 on GitHub](https://github.com/DerrickWood/kraken2/)  |
        | Software Documentation | <https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown> |
        | Original Publication(s) | [Improved metagenomic analysis with Kraken 2](https://link.springer.com/article/10.1186/s13059-019-1891-0) |
