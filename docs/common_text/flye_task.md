??? task "`flye`"

    Flye is a _de novo_ assembler for long read data using repeat graphs. Compared to de Bruijn graphs, which require exact k-mer matches, repeat graphs can use approximate matches which better tolerates the error rate of ONT data.

<!-- if: theiaviral -->
    It can be enabled by setting the `call_raven` parameter to `false`. The `flye` task is used as a fallback option if the `raven` task fails during execution (see task `raven` for more details).
<!-- endif -->

    ??? dna "`flye_read_type`" 
        This input parameter specifies the type of sequencing reads being used for assembly. This parameter significantly impacts the assembly process and should match the characteristics of your input data. Below are the available options:
    
        | **Parameter** | **Explanation** |
        | --- | --- |
        | `--nano-hq` (default) | Optimized for ONT high-quality reads, such as Guppy5+ SUP or Q20 (<5% error). Recommended for ONT reads processed with Guppy5 or newer |
        | `--nano-raw` | For ONT regular reads, pre-Guppy5 (<20% error) |
        | `--nano-corr` | ONT reads corrected with other methods (<3% error) |
        | `--pacbio-raw` | PacBio regular CLR reads (<20% error) |
        | `--pacbio-corr` | PacBio reads corrected with other methods (<3% error) |
        | `--pacbio-hifi` | PacBio HiFi reads (<1% error) |
    
        Refer to the Flye documentation for detailed guidance on selecting the appropriate `flye_read_type` based on your sequencing data and additional optional paramaters.

<!-- if: theiaviral -->
    ???+ warning "Important"
        In this workflow, de novo assembly is used solely to facilitate the selection of a closely related reference genome. If the user provides an input `reference_fasta`, all subsequent assembly and reference selections tasks will be skipped, including:

        - `flye`
        - `checkv_denovo`
        - `quast_denovo`
        - `skani`
        - `ncbi_datasets`
<!-- endif -->

    ???+ warning "Non-deterministic output(s)"
        This task may yield non-deterministic outputs.

    !!! techdetails "Flye Technical Details"
        |  | Links |
        | --- | --- |
        | WDL Task | [task_flye.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_flye.wdl) |
        | Software Source Code | [Flye on GitHub](https://github.com/fenderglass/Flye) |
        | Software Documentation | [Flye Documentation](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md) |
        | Original Publication(s) | [Assembly of long, error-prone reads using repeat graphs](https://www.nature.com/articles/s41587-019-0072-8) |