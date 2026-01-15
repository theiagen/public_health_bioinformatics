??? task "`TBProfiler`: Lineage and Drug Susceptibility Prediction ==_for Illumina and ONT only_=="
    [TBProfiler](https://github.com/jodyphelan/TBProfiler) identifies _Mycobacterium tuberculosis_ complex species, lineages, sub-lineages and drug resistance-associated mutations.

    TBProfiler aligns the input reads to the H37Rv (NC_000962.3/AL123456.3) reference genome with `bwa mem` (or `minimap2`) and then calls variants using `gatk` (default), though other options are available (`bcftools`, `freebayes`, `lofreq`, `pilon`). After mutations are called and filtered, they are compared against [TBProfiler's database (TBDB)](https://github.com/jodyphelan/tbdb).

    A number of outputs are made available from TBProfiler, all of which can be found in the summary results JSON file. Although the JSON file contains the most information, it is not very human readable, which is why the results CSV and TXT files have been made available to the user as outputs. TBProfiler is able to detect the identified lineage and any sublineages, determine the predicted type of drug resistance, a lists the genes that are associated with resistance, along with several other useful outputs.     

    !!! techdetails "TBProfiler Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_tbprofiler.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/mycobacterium/task_tbprofiler.wdl) |
        | Software Source Code | [TBProfiler on GitHub](https://github.com/jodyphelan/TBProfiler) |
        | Software Documentation | [TBProfiler Docs on GitHub](https://jodyphelan.github.io/tb-profiler-docs/en/) |
        | Original Publication(s) | [Integrating informatics tools and portable sequencing technology for rapid detection of resistance to anti-tuberculous drugs](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0650-x) |
