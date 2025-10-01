??? task "`Bakta`: Assembly Annotation (alternative)"

    To activate this task, set `genome_annotation` to `"bakta"`.

    `Bakta` is intended for the annotation of bacterial genomes and plasmids and is used to identify and describe regions of interest within the genome.

    In addition to the standard annotation outputs, `Bakta` also provides a plot summarizing the annotation results, which can be useful for visualizing genome features.

    **Bakta Database Options**

    Our implementation of `Bakta` supports three database configurations:

    - **Light**: Optimized for faster performance and lower resource usage, with a focused set of core reference data for most bacterial genome annotations. Recommended for quick annotations or limited computational resources. Specify "light" for the `bakta_db` input.
    - **Full** (default): Comprehensive with extensive reference annotations, suitable for detailed and accurate annotations. Specify "full" for the `bakta_db` input.
    - **Custom**: Allows users to provide a Bakta-compatible database stored in Google Cloud Storage Bucket. This file must be a .tar.gz archive containing a properly formatted Bakta database with a valid version.json. Please see the [Bakta database documentation](https://github.com/oschwengers/bakta#database) for detailed formatting requirements. Example: `"bakta_db": "gs://my-bucket/custom_bakta_db.tar.gz"`

    !!! techdetails "Bakta Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_bakta.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/annotation/task_bakta.wdl) |
        | Software Source Code | [Bakta on GitHub](https://github.com/oschwengers/bakta) |
        | Software Documentation | [Bakta on GitHub](https://github.com/oschwengers/bakta) |
        | Original Publication(s) | [Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000685) |
