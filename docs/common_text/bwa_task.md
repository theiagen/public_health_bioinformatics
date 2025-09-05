<!-- if: theiameta -->
??? task "`bwa`: Read alignment to the assembly"
    ==If a reference is _not_ provided==, BWA (Burrow-Wheeler Aligner) is used to align the clean reads to the Pilon-polished assembly_fasta.
<!-- endif -->
<!-- if: digger -->
??? task "`bwa`: Read Alignment to the Assembly"
    BWA (Burrow-Wheeler Aligner) is used to align the cleaned read files to generated assembly file in order to generate an alignment. The resulting BAM file is directly passed to the Pilon task to polish the assembly for errors.

<!-- if: freyja -->
??? task "`bwa` Details"
    This task aligns the cleaned short reads (Illumina) to the reference genome provided by the user.
<!-- endif -->
<!-- if: theiaviral -->
??? task "`bwa`"

    The `bwa` task is a wrapper for the BWA alignment tool. It utilizes the BWA-MEM algorithm to map cleaned reads to the reference genome, either selected by the `skani` task or provided by the user input `reference_fasta`. This creates a BAM file which is then sorted using the command `samtools sort`.
<!-- endif -->
    !!! techdetails "BWA Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_bwa.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/alignment/task_bwa.wdl) |
        | Software Source Code | [BWA on GitHub](https://github.com/lh3/bwa) |
        | Software Documentation | [BWA Documentation](https://bio-bwa.sourceforge.net/) |
        | Original Publication(s) | [Fast and accurate short read alignment with Burrows-Wheeler transform](https://doi.org/10.1093/bioinformatics/btp324) |
