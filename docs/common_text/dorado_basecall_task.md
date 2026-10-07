---
title: Task Fragment `dorado_basecall`
fragment: true
---
??? task "`dorado_basecall`: Basecalling POD5 files"
    The basecalling task takes POD5 files as input and converts each individual POD5 file into 'BAM' format using either the default or user-specified model. This step leverages GPU acceleration for efficient processing.

    Please see the [Dorado documentation](https://dorado-docs.readthedocs.io/en/latest/basecaller/basecall_overview/) for more details, but what follows is a brief overview of the basecalling process:

    1. POD5 files are pre-processed via signal scaling and normalization.
    2. The machine learning algorithm decodes the sequence signals into nucleotide base calls. There are different machine learning models that can be specified as input; see [Model Type Selection](#model-type-selection).
    3. [Barcode classification](https://dorado-docs.readthedocs.io/en/latest/barcoding/barcoding/) is performed based on the indicated kit name to enable downstream demultiplexing.

        !!! info "Barcode Trimming"
            Barcode trimming is purposefully **disabled** during the basecalling step to ensure accurate demultiplexing in subsequent workflow steps.

    4. Modified basecalling can be performed if indicated through [modification to the model name](https://dorado-docs.readthedocs.io/en/latest/basecaller/mods/).
    5. [Reads are split](https://dorado-docs.readthedocs.io/en/latest/basecaller/read_splitting/) when a single read contains multiple concatenated reads.

    Other options are available, but currently Dorado_Basecalling_PHB does not support them. Please contact <support@theiagen.com> if you would like additional options.

    !!! techdetails "Dorado Basecalling Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_dorado_basecall.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/basecalling/task_dorado_basecall.wdl) |
        | Software Source Code | [Dorado on GitHub](https://github.com/nanoporetech/dorado/) |
        | Software Documentation | [Dorado ReadTheDocs](https://dorado-docs.readthedocs.io/en/latest/) |
