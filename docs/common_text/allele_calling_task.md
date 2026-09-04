---
title: Task Fragment `allele_calling`
fragment: true
---

??? task "`allele_calling`: PulseNet 2.0 Hash-Based Allele Calling"
    The Allele Calling module is used by PulseNet 2.0 to generate MLST calls for core and accessory genes with corresponding quality metrics for certain supported organisms.

    The module works by comparing genome assemblies against [reference allele sequences](https://github.com/ncezid-biome/pn2.0-mlst-databases/tree/main) using a BLASTn approach to find the presence of each locus. The query allele sequence is defined by the presence of start and stop codons (without nonsense mutations) and a similarity threshold (organism-specific) against a reference allele. Loci that were likely repeated (fully or partially) elsewhere in the genome are ignored. The query sequences are hashed using the 64-bit MD5 algorithm and then transformed into a 56-bit integer.

    These results are stored into multiple output files, for core and accessory MLST calls, alongside quality metrics such as core genome percentage and a pass/fail flag.

    The following default options are used:

<!-- if: campy -->
    - similarity = 70
    - database = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/CAMPY.tar.gz"
<!-- endif -->
<!-- if: cbot -->
    - similarity = 85
    - database = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/CBOT.tar.gz"
<!-- endif -->
<!-- if: crono -->
    - similarity = 80
    - database = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/CRONO.tar.gz"
<!-- endif -->
<!-- if: stec -->
    - similarity = 85
    - database = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/STEC.tar.gz"
<!-- endif -->
<!-- if: listeria -->
    - similarity = 85
    - database = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/LISTERIA.tar.gz"
<!-- endif -->
<!-- if: salm -->
    - similarity = 75
    - database = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/SALM.tar.gz"
<!-- endif -->
<!-- if: vibrio -->
    - similarity = 85
    - database = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/VIBR.tar.gz"
<!-- endif -->
<!-- if: yersinia -->
    - similarity = 85
    - database = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/YERSINIA.tar.gz"
<!-- endif -->

    These parameters, and others, are set in a preceeding [`allele_calling_parameters` subworkflow](https://github.com/theiagen/public_health_bioinformations/blob/main/workflows/utilities/wf_allele_calling_parameters.wdl). These parameters include the organism-specific databases sourced from the [pn2.0-mlst-databases](https://github.com/ncezid-biome/pn2.0-mlst-databases/tree/main) repository, specific pathing information unique to each database, and [organism-specific similarity thresholds](https://github.com/ncezid-biome/pulsenet2.0-bfx/blob/47644186f2df27e9f01a000d47c451135a75f65d/main/utils/utils.nf#L40).

    !!! techdetails "Allele Calling Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_allele_calling.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/multi/task_allele_calling.wdl) |
        | Software Source Code | [PulseNet 2.0 BFX AlleleCalling Process](https://github.com/ncezid-biome/pulsenet2.0-bfx/tree/main/processes/AlleleCalling) |
        | Software Documentation |  [PulseNet 2.0 BFX AlleleCalling Process on GitHub](https://github.com/ncezid-biome/pulsenet2.0-bfx/tree/main/processes/AlleleCalling#readme) |
