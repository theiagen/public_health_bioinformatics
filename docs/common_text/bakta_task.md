??? task "`Bakta`: Assembly Annotation (alternative)"

    Assembly annotation is available via `Bakta` as an alternative to `Prokka`. When `Bakta` annotation is used, `Prokka` is not.

    `Bakta` is intended for annotation of Bacteria and plasmids only, and is best described [here](https://github.com/oschwengers/bakta#description)!

    In addition to the standard annotation outputs, `Bakta` also provides a plot summarizing the annotation results, which can be useful for visualizing genome features.

    **Bakta Database Options**

    `Bakta` supports three database configurations:

    **Light** Database: Optimized for faster performance and lower resource usage, with a focused set of core reference data for most bacterial genome annotations. Recommended for quick annotations or limited computational resources. Specify "light" for the `bakta_db` input.

    **Full** Database (default): Comprehensive with extensive reference annotations, suitable for detailed and accurate annotations. Specify "full" for the `bakta_db` input.

    **Custom** Database: Allows users to provide a Bakta-compatible database stored in Google Cloud Storage Must be a .tar.gz archive containing a properly formatted Bakta database with a valid version.json Follow the [Bakta database documentation](https://github.com/oschwengers/bakta#database) for detailed formatting requirements. Example: `"bakta_db": "gs://my-bucket/custom_bakta_db.tar.gz"`

    `Bakta` also supports the use of user-provided proteins via `proteins` input. These are used in the functional annotation process. See the [User-provided protein sequences](https://github.com/oschwengers/bakta/blob/main/README.md#user-provided-protein-sequences) documentation for further details. 

    !!! techdetails "Bakta Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_bakta.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/annotation/task_bakta.wdl) |
        | Software Source Code | [bakta](https://github.com/oschwengers/bakta) |
        | Software Documentation | <https://github.com/oschwengers/bakta> |
        | Original Publication(s) | [Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000685) |
