??? task "`ivar consensus`"

    The `consensus` task wraps the [iVar](https://andersen-lab.github.io/ivar/html/index.html) tool to generate a reference-based consensus assembly from the sorted BAM file produced by the `bwa` task. It uses the `ivar consensus` command to call variants and generate a consensus sequence based on those mapped reads. The `consensus` task will filter all variant calls based on user-defined parameters, including `min_map_quality`, `min_depth`, and `min_allele_freq`. This task will return a consensus sequence in FASTA format and the samtools mpileup output.

    This task is functional for segmented viruses by iteratively executing iVar on a contig-by-contig basis and concantenating resulting consensus contigs.

    ??? dna "`min_depth` input parameter"
        This parameter accepts an integer value to set the minimum read depth for variant calling and subsequent consensus sequence generation. The default value is `10`.

    ??? dna "`min_map_quality` input parameter"
        This parameter accepts an integer value to set the minimum mapping quality for variant calling and subsequent consensus sequence generation. The default value is `20`.

    ??? dna "`min_allele_freq` input parameter"
        This parameter accepts a float value to set the minimum allele frequency for variant calling and subsequent consensus sequence generation. The default value is `0.6`.

    !!! techdetails "iVar Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_ivar_consensus.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_ivar_consensus.wdl) |
        | Software Source Code | [Ivar on GitHub](https://andersen-lab.github.io/ivar/html/) |
        | Software Documentation | [Ivar Documentation](https://andersen-lab.github.io/ivar/html/manualpage.html) |
        | Original Publication(s) | [An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar](http://dx.doi.org/10.1186/s13059-018-1618-7) |