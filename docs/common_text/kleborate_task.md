??? task "`Kleborate`: Species Identification, MLST, Serotyping, AMR and Virulence Characterization"
    [Kleborate](https://github.com/katholt/Kleborate) is a tool to identify the _Klebsiella_ species, MLST sequence type, serotype, virulence factors (ICE_Kp_ and plasmid associated), and AMR genes and mutations. Serotyping is based on the capsular (K antigen) and lipopolysaccharide (LPS) (O antigen) genes. The acquired resistance genes identified by Kleborate can be found in [the Kleborate documentation here](https://kleborate.readthedocs.io/en/latest/kpsc_modules.html#acquired-amr-genes), along with other useful information regarding all of Kleborate's modules.

    `Kaptive` can be run as well by setting the optional input variable `kleborate_skip_kaptive` to `false` in order to recieve the K antigen and O antigen locus typing via _wzi_ alleles.

    !!! techdetails "Kleborate Technical Details"        
        |  | Links |
        | --- | --- |
        | Task | [task_kleborate.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/klebsiella/task_kleborate.wdl) |
        | Software Source Code | [Kleborate on GitHub](https://github.com/klebgenomics/Kleborate) |
        | Software Documentation | [Kleborate Documentation on ReadTheDocs](https://kleborate.readthedocs.io/en/latest/) |
        | Original publication | [A genomic surveillance framework and genotyping tool for Klebsiella pneumoniae and its related species complex](https://www.nature.com/articles/s41467-021-24448-3)<br>[Identification of Klebsiella capsule synthesis loci from whole genome data](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000102)<br> |
