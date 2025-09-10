??? task "`ivar_consensus`: Consensus Assembly"
    iVar's `consensus` tool generates a reference-based consensus assembly. Several parameters can be set that determine the stringency of the consensus assembly, including minimum quality, minimum allele frequency, and minimum depth.

<!-- if: theiacov -->
    For TheiaCoV, the following default parameters are used:
    
    - minimum quality: 20
    - minimum depth: 100
    - minimum allele frequency: 0.6
<!-- endif -->

<!-- if: theiaviral -->
    This task is functional for segmented viruses by iteratively executing iVar on a contig-by-contig basis and concantenating resulting consensus contigs.

    ??? dna "`min_depth` input parameter"
        This parameter accepts an integer value to set the minimum read depth for variant calling and subsequent consensus sequence generation. The default value is `10`.

    ??? dna "`min_map_quality` input parameter"
        This parameter accepts an integer value to set the minimum mapping quality for variant calling and subsequent consensus sequence generation. The default value is `20`.

    ??? dna "`min_allele_freq` input parameter"
        This parameter accepts a float value to set the minimum allele frequency for variant calling and subsequent consensus sequence generation. The default value is `0.6`.
<!-- endif -->

    !!! techdetails "iVar Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_ivar_consensus.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_ivar_consensus.wdl) |
        | Software Source Code | [Ivar on GitHub](https://andersen-lab.github.io/ivar/html/) |
        | Software Documentation | [Ivar Documentation](https://andersen-lab.github.io/ivar/html/manualpage.html) |
        | Original Publication(s) | [An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar](http://dx.doi.org/10.1186/s13059-018-1618-7) |
