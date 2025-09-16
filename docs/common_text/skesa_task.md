??? task "`SKESA`: _De novo_ Assembly (default)"
    This task is activated by default.

    `SKESA` (Strategic K-mer Extension for Scrupulous Assemblies) is a _de novo_ assembler that is fairly conservative and introduces breaks in the genome at repeat regions. This leads to higher sequence quality but more fragmented assemblies, which, depending on the final analysis goal, can be either highly preferred or detrimental. Designed for Illumina reads and haploid genomes, SKESA is the default assembler in the `digger_denovo` subworkflow.

    !!! techdetails "SKESA Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_skesa.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_skesa.wdl) |
        | Software Source Code | [SKESA on GitHub](https://github.com/ncbi/SKESA) |
        | Software Documentation | [SKESA on GitHub](https://github.com/ncbi/SKESA) |
        | Original Publication(s) | [SKESA: strategic k-mer externsion for scrupulous assemblies](https://doi.org/10.1186/s13059-018-1540-z) |
