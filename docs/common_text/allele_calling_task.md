---
title: Task Fragment `allele_calling`
fragment: true
---

??? task "`allele_calling`: PulseNet 2.0 Hash-Based Allele Calling"
    The Allele Calling module is used by PulseNet 2.0 to generate MLST calls for core and accessory genes with corresponding quality metrics for certain supported organisms.

    The module works by comparing genome assemblies against [reference allele sequences](https://github.com/ncezid-biome/pn2.0-mlst-databases/tree/main) using a BLASTn approach to find the presence of each locus. The query allele sequence is defined by the presence of start and stop codons (without nonsense mutations) and a similarity threshold (organism-specific) against a reference allele. Loci that were likely repeated (fully or partially) elsewhere in the genome are ignored. The query sequences are hashed using the 64-bit MD5 algorithm and then transformed into a 56-bit integer.

    These results are stored into multiple output files, for core and accessory MLST calls, alongside quality metrics such as core genome percentage and a pass/fail flag.

    Organism-specific parameters are set in a preceeding `allele_calling_parameters` subworkflow. These parameters include the organism-specific databases sourced from the [pn2.0-mlst-databases](https://github.com/ncezid-biome/pn2.0-mlst-databases/tree/main) repository, specific pathing information unique to each database, and [organism-specific similarity thresholds](https://github.com/ncezid-biome/pulsenet2.0-bfx/blob/47644186f2df27e9f01a000d47c451135a75f65d/main/utils/utils.nf#L40).

    !!! techdetails "Allele Calling Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_allele_calling.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/multi/task_allele_calling.wdl) |
        | Software Source Code | [PulseNet 2.0 BFX AlleleCalling Process](https://github.com/ncezid-biome/pulsenet2.0-bfx/tree/main/processes/AlleleCalling) |
        | Software Documentation |  [PulseNet 2.0 BFX AlleleCalling Process on GitHub](https://github.com/ncezid-biome/pulsenet2.0-bfx/tree/main/processes/AlleleCalling#readme) |
