??? task "`SISTR`: _Salmonella_ Serovar Prediction"
    [SISTR](https://github.com/phac-nml/sistr_cmd) performs _Salmonella spp_. serotype prediction using antigen gene and cgMLST gene alleles. In TheiaProk, SISTR is run on genome assemblies and uses the default database setting (which contains smaller "centroid" alleles or representative alleles instead of the full set of cgMLST alleles). The full set of cgMLST alleles can be activated by setting the `sistr_use_full_cgmlst_db` optional input variable to `true`. It also runs a QC module to determine the level of confidence in the serovar prediction (please see [the section on QC in the SISTR documentation here](https://github.com/phac-nml/sistr_cmd#qc-by-sistr_cmd---qc)).

    SISTR uses a database of _Salmonella_ serovar determination antigens, cgMLST profiles, and MASH sketches of appropriate reference genomes. BLAST is used to compare input assemblies to this database for serotyping, and the [Mash MinHash](https://mash.readthedocs.io/en/latest/) algorithm is used for serovar prediction.

    !!! techdetails "SISTR Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_sistr.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/salmonella/task_sistr.wdl) |
        | Software Source Code | [SISTR on GitHub](https://github.com/phac-nml/sistr_cmd) |
        | Software Documentation | [SISTR on GitHub](https://github.com/phac-nml/sistr_cmd) |
        | Original Publication(s) | [The _Salmonella In Silico_ Typing Resource (SISTR): an open web-accessible tool for rapidly typing and subtyping draft _Salmonella_ genome assemblies.](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0147101) |
