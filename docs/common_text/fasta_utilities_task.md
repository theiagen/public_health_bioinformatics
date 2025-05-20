??? task "`fasta_utilities`"

    The `fasta_utilities` task utilizes samtools to index a reference fasta file.
<!-- if: theiaviral-->
    This reference is selected by the `skani` task or provided by the user input `reference_fasta`. This indexed reference genome is used for downstream variant calling and consensus generation tasks.
<!-- endif -->

    !!! techdetails "`fasta_utilities` Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_fasta_utilities.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/data_handling/task_fasta_utilities.wdl) |
        | Software Source Code | [samtools on GitHub](https://github.com/samtools/samtools) |
        | Software Documentation | [samtools](https://www.htslib.org/doc/samtools.html) |
        | Original Publication(s) | [The Sequence Alignment/Map format and SAMtools](https://doi.org/10.1093/bioinformatics/btp352)<br>[Twelve Years of SAMtools and BCFtools](https://doi.org/10.1093/gigascience/giab008) |