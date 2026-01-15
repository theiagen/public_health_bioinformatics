??? task "`raven`"

    The `raven` task is used to create a *de novo* assembly from cleaned reads. Raven is an overlap-layout-consensus based assembler that accelerates the overlap step, constructs an assembly graph from reads pre-processed with pile-o-grams, applies a novel and robust graph simplification method based on graph drawings, and polishes unambiguous graph paths using Racon. 

<!-- if: theiaviral -->
    Based on internal benchmarking against Flye and results reported by [Cook et al. (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11092197/), Raven is faster, produces more contiguous assemblies, and yields more complete genomes within TheiaViral according to CheckV quality assessment (see task `checkv` for technical details).

    ??? dna "`call_raven` input parameter"
        This parameter controls whether or not the `raven` task is called by the workflow. By default, `call_raven` is set to `true` because Raven is used as the primary assembler. Raven is generally recommended for most users, but it might not perform optimally on all datasets. If users encounter issues with Raven, they can set the `call_raven` variable to `false` to bypass the `raven` task and instead *de novo* assemble using Flye (see task `flye` for details). Additionally, if the Raven task fails during execution, the workflow will automatically fall back to using Flye for *de novo* assembly.
<!-- endif -->

    ???+ warning "Error traceback"
        Raven may fail with cryptic "segmentation fault" (segfault) errors or by failing to output an output file. It is difficult to traceback the source of these issues, though increasing the `memory` parameter may resolve some errors.

    ???+ warning "Non-deterministic output(s)"
        This task may yield non-deterministic outputs.

    !!! techdetails "Raven Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_raven.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_raven.wdl) |
        | Software Source Code | [Raven on GitHub](https://github.com/lbcb-sci/raven)
        | Software Documentation | [Raven Documentation](https://github.com/lbcb-sci/raven/blob/master/README.md)
        | Original Publication(s) | [Time- and memory-efficient genome assembly with Raven](https://doi.org/10.1038/s43588-021-00073-4) |