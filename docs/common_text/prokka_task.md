??? task "`Prokka`: Assembly Annotation (default)"

    Assembly annotation is available via `Prokka` as default.

    [`Prokka`](https://github.com/tseemann/prokka) is a prokaryotic genome annotation tool used to identify and describe features of interest within the genome sequence. Prokka annotates the genome by querying three core databases: [ISfinder](https://isfinder.biotoul.fr/), [NCBI's Bacterial Antimicrobial Resistance Reference Gene Database](https://www.ncbi.nlm.nih.gov/bioproject/313047), and [UniProtKB (SwissProt)](https://www.uniprot.org/uniprotkb?query=reviewed%3Atrue). Additional databases can be used or specified, with instructions on how to do so [located in the Prokka README](https://github.com/tseemann/prokka#databases).

    The most versatile output from Prokka is likely the GFF3 file, as it contains all of the generated information, though other file formats are available for the instances when a reduction of the GFF3 is useful.

    !!! techdetails "Prokka Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_prokka.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/annotation/task_prokka.wdl) |
        | Software Source Code | [Prokka on GitHub](https://github.com/tseemann/prokka) |
        | Software Documentation | [Prokka on GitHub](https://github.com/tseemann/prokka) |
        | Original Publication(s) | [Prokka: rapid prokaryotic genome annotation](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517?login=false) |
