---
title: Task Fragment `bwa`
fragment: true
---
??? task "`BWA`: Read Alignment to the Assembly"
    
<!-- if: theiameta -->
    ==If a reference is _not_ provided==, BWA (Burrow-Wheeler Aligner) is used to align the cleaned reads to the Pilon-polished assembly_fasta.
<!-- endif -->
<!-- if: digger -->
    BWA (Burrow-Wheeler Aligner) is used to align the cleaned read files to a generated assembly file. The resulting BAM file is directly passed to the Pilon task to polish the assembly for errors.
<!-- endif -->

<!-- if: theiacov -->
    BWA (Burrow-Wheeler Aligner) is used to align the cleaned read files to a reference genome, either determined by the user or provided by the [_organism-specific parameters_ section](./theiacov.md#org-specific) (see above). The resulting BAM file is used for primer trimming, variant calling, and consensus generation in downstream tasks.
<!-- endif -->

<!-- if: freyja -->
    BWA (Burrow-Wheeler Aligner) is used to align the cleaned read files to a reference genome provided by the user.
<!-- endif -->
<!-- if: theiaviral -->
    BWA (Burrow-Wheeler Aligner) is used to align the cleaned read files to a reference genome either selected by skani or provided by the user with the `reference_fasta` input. This creates a BAM file which is then sorted using samtools.
<!-- endif -->
    !!! techdetails "BWA Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_bwa.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/alignment/task_bwa.wdl) |
        | Software Source Code | [BWA on GitHub](https://github.com/lh3/bwa) |
        | Software Documentation | [BWA Documentation](https://bio-bwa.sourceforge.net/) |
        | Original Publication(s) | [Fast and accurate short read alignment with Burrows-Wheeler transform](https://doi.org/10.1093/bioinformatics/btp324) |
