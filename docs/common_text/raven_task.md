??? task "`raven`"

    The `raven` task is used to create a de novo assembly from cleaned reads. Raven is an overlap-layout-consensus based assembler that accelerates the overlap step, constructs an assembly graph from reads pre-processed with pile-o-grams, applies a novel and robust graph simplification method based on graph drawings, and polishes unambiguous graph paths using Racon. 

<!-- if: theiaviral -->
    Based on comprehensive benchmarking against Flye and results reported by [Cook et al. (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11092197/), Raven is faster, produces more contiguous assemblies, and yields more complete genomes according to CheckV metrics (see task `checkv` for technical details).

    ??? dna "`skip_raven`"
        This parameter controls whether or not the `raven` task is skipped by the workflow. By default, `skip_raven` is set to `false` because Raven is used as the primary assembler. Raven is generally recommended for most users, but it might not perform optimally on all datasets. If users encounter issues with Raven, they can set the `skip_raven` variable to `true` to bypass the `raven` task and instead de novo assemble using Flye (see task `flye` for details). Additionally, if the raven task ever fails during execution, the workflow will automatically fall back to using Flye for de novo assembly.

    ???+ warning "Important"
        In this workflow, de novo assembly is used solely to facilitate the selection of a closely related reference genome. If the user provides an input `reference_fasta`, all subsequent assembly and reference selections tasks will be skipped, including:

        - `raven`
        - `flye`
        - `checkv_denovo`
        - `quast_denovo`
        - `skani`
        - `ncbi_datasets`
<!-- endif -->

    !!! techdetails "Raven Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_raven.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_raven.wdl) |
        | Software Source Code | [Raven on GitHub](https://github.com/lbcb-sci/raven)
        | Software Documentation | [Raven Documentation](https://github.com/lbcb-sci/raven/blob/master/README.md)
        | Original Publication(s) | [Raven Paper](https://doi.org/10.1038/s43588-021-00073-4) |