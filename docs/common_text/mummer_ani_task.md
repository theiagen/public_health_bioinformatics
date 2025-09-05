
??? task "`MUMmer_ANI`: Taxon Assignment using Average Nucleotide Identity (optional)"
    To activate this task, set `call_ani` to `true`.

    Average Nucleotide Identity (ANI) is a useful approach for taxonomic identification. The higher the percentage ANI of a query sequence to a given reference genome, the more likely the sequence is the same taxa as the reference. 

    ANI is calculated in TheiaProk using [a perl script written by Lee Katz](https://github.com/lskatz/ani-m) (ani-m.pl). This uses [MUMmer](http://mummer.sourceforge.net/) to rapidly align entire query assemblies to one or more reference genomes. By default, TheiaProk uses a set of 43 reference genomes in [RGDv2](https://github.com/StaPH-B/docker-builds/blob/master/build-files/fastani/1.34-RGDV2/RGDv2-metadata.tsv), a database containing genomes of enteric pathogens commonly sequenced by CDC EDLB & PulseNet participating laboratories. The user may also provide their own reference genome. After genome alignment with MUMmer, ani-m.pl calculates the average nucleotide identity and percent bases aligned between 2 genomes (query and reference genomes)

    The default database of reference genomes used is called "Reference Genome Database version 2" AKA "RGDv2". This database is composed of 43 enteric bacteria representing 32 species and is intended for identification of enteric pathogens and common contaminants. It contains six Campylobacter spp., three Escherichia/Shigella spp., one *Grimontia hollisae*, six *Listeria spp.*, one *Photobacterium damselae*, two *Salmonella spp.*, and thirteen *Vibrio spp.* 

    2 Thresholds are utilized to prevent false positive hits. The `ani_top_species_match` will only report a genus & species match if both thresholds are surpassed. Both of these thresholds are set to match those used in BioNumerics for PulseNet organisms.

    1. `ani_threshold` default value of 80.0
    2. `percent_bases_aligned_threshold` default value of 70.0

    For more information on RGDv2 database of reference genomes, please see [the publication here](https://doi.org/10.3389/fmicb.2023.1225207).

    !!! techdetails "MUMmer_ANI Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_mummer_ani.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/advanced_metrics/task_mummer_ani.wdl) |
        | Software Source Code | [ani-m on GitHub](https://github.com/lskatz/ani-m)<br>[MUMmer on GitHub](https://github.com/mummer4/mummer) |
        | Software Documentation | [ani-m on GitHub](https://github.com/lskatz/ani-m)<br>[MUMmer on SourceForge](https://mummer.sourceforge.net/) |
        | Original Publication(s) | _MUMmer4_: [MUMmer4: A fast and versatile genome alignment system](https://doi.org/10.1371/journal.pcbi.1005944)<br>_RGDv2 database_: [Rapid identification of enteric bacteria from whole genome sequences using average nucleotide identity metrics](https://doi.org/10.3389/fmicb.2023.1225207) |
