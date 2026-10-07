---
title: Task Fragment `dorado_demux`
fragment: true
---
??? task "`dorado_demux`: Produces barcode-specific FASTQ files"
    This task takes each basecalled BAM file and demultiplexes it based on the identified barcodes found during basecalling. An individual FASTQ file is generated for each barcode found per BAM file. All FASTQ files that are associated with a single barcode are then merged.

    !!! info "Disabling Barcode Trimming"
        By default, barcodes _are_ trimmed during demultiplexing.

        This can be disabled by setting the optional input variable `demux_no_trim` to `true`. This allows users to retain untrimmed reads for troubleshooting, such as inspecting reads in the "unclassified" folder when reads are mis-binned or other data issues occur.

    !!! techdetails "Dorado Demultiplexing Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_dorado_demux.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/basecalling/task_dorado_demux.wdl) |
        | Software Source Code | [Dorado on GitHub](https://github.com/nanoporetech/dorado/) |
        | Software Documentation | [Dorado ReadTheDocs](https://dorado-docs.readthedocs.io/en/latest/) |
