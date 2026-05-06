---
title: Task Fragment `fasta_utilities`
fragment: true
---
??? task "`fasta_utilities`: Reference Indexing"
    The `fasta_utilities` task uses SAMtools to index a reference fasta file.
    
<!-- if: theiaviral-->
    This reference is selected by the `skani` task or provided by the user input `reference_fasta`. This indexed reference genome is used for downstream variant calling and consensus generation tasks.
<!-- endif -->

    !!! techdetails "`fasta_utilities` Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_fasta_utilities.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/data_handling/task_fasta_utilities.wdl) |
        | Software Source Code | [SAMtools on GitHub](https://github.com/samtools/samtools) |
        | Software Documentation | [SAMtools](https://www.htslib.org/doc/samtools.html) |
        | Original Publication(s) | [The Sequence Alignment/Map format and SAMtools](https://doi.org/10.1093/bioinformatics/btp352)<br>[Twelve Years of SAMtools and BCFtools](https://doi.org/10.1093/gigascience/giab008) |
