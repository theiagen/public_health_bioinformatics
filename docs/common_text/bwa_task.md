<!-- if: theiameta -->
??? task "`bwa`: Read alignment to the assembly"
    ==If a reference is _not_ provided==, BWA (Burrow-Wheeler Aligner) is used to align the clean reads to the Pilon-polished assembly_fasta.
<!-- endif -->
<!-- if: freyja -->
??? task "`bwa` Details"
    This task aligns the cleaned short reads (Illumina) to the reference genome provided by the user.
<!-- endif -->

    !!! techdetails "BWA Technical Details"
    
        |  | Links |
        | --- | --- |
        | Task | [task_bwa.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/alignment/task_bwa.wdl) |
        | Software Source Code | <https://github.com/lh3/bwa> |
        | Software Documentation | <https://bio-bwa.sourceforge.net/> |
        | Original Publication(s) | [Fast and accurate short read alignment with Burrows-Wheeler transform](https://doi.org/10.1093/bioinformatics/btp324) |