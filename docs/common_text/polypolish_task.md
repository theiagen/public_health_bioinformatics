??? task "`Polypolish`: Hybrid Assembly Polishing ==_for ONT and Illumina data_=="
    If short reads are provided with the optional `illumina_read1` and `illumina_read2` inputs, Polypolish will use those short-reads to correct errors in the long-read assemblies. Uniquely, Polypolish uses the short-read alignments where each read is aligned to _all_ possible locations, meaning that even repeat regions will have error correction.

    !!! techdetails "Polypolish Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_polypolish.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/polishing/task_polypolish.wdl) |
        | Software Source Code | [Polypolish on GitHub](https://github.com/rrwick/Polypolish) |
        | Software Documentation | [Polypolish Documentation](https://github.com/rrwick/Polypolish#polypolish) |
        | Original Publication(s) | [Polypolish: short-read polishing of long-read bacterial genome assemblies](https://doi.org/10.1371/journal.pcbi.1009802)<br>[How low can you go? Short-read polishing of Oxford Nanopore bacterial genome assemblies](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001254) |
