??? task "`BBDuk`: Adapter Trimming and PhiX Removal"
    Adapters are manufactured oligonucleotide sequences attached to DNA fragments during the library preparation process. In Illumina sequencing, these adapter sequences are required for attaching reads to flow cells. You can read more about [Illumina adapters here](https://emea.support.illumina.com/bulletins/2020/06/illumina-adapter-portfolio.html). For genome analysis, it's important to remove these sequences since they're not actually from your sample. If you don't remove them, the downstream analysis may be affected.

    By default, the `bbduk` task will:

    - **Repair disordered read pairs** (if they exist) so that the first read in file 1 is the same mate of the first read in file2. See the [`Repair Guide`](https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/repair-guide/) from the [BBTools](https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/) package.

    - **Remove PhiX contamination** by filtering out all reads that have a 31-mer match to [PhiX](https://emea.illumina.com/products/by-type/sequencing-kits/cluster-gen-sequencing-reagents/phix-control-v3.html). PhiX is a viral genome that is often used as a control in Illumina sequencing runs. Removing PhiX sequences helps to ensure that the data reflects only the target organism's genome. By default this task uses the built-in PhiX reference fasta provided with BBTools (see [here](https://github.com/bbushnell/BBTools/blob/master/resources/phix174_ill.ref.fa.gz)), but a custom PhiX reference can be provided via the `phix_fasta` input parameter.

    - **Trim adapters** with [BBDuk](https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/)  (*"Bestus Bioinformaticus" Decontamination Using Kmers*). By default this uses the built-in adapter sequences provided with BBTools (see [here](https://github.com/bbushnell/BBTools/blob/master/resources/adapters.fa)). If you have custom adapter sequences specific to your library preparation, you can provide them via the `adapters_fasta` input parameter.

<!-- if: theiaviral -->
    Optionally, primers can be trimmed using an alignment-based approach implemented in BBDuk. To activate this functionality, primers must be provided via a standard format FASTA file (`primers_fasta`), with each entry depicting a separate primer, or can be inputted as a comma-delimited string of primers (`primers_literal`).

    By default, the kmer size used for primer trimming will automatically be set based on the individual lengths of each primer sequence provided. Additionally primers will be trimmed from both the 5' and 3' ends of reads.

    The sensitivity of alignment-based primer trimming can be tuned via several parameters:

    - `primers_hamming_distance`: This variable indicates the number of allowed mismatched bases during the alignment of primers to reads. Default is `1`.
    - `primers_max_start_offset`: This allows the trimming of primers to start "X" number of bases into the read. Default is `10`.
    - `primers_mask_middle`: This is a boolean input that, when set to "true", will mask the middle base of a primer sequence during alignment, which can be useful for handling degenerate primers. Default is `"false"`.
    - `primers_reverse_complement` This is a boolean input that, when set to "true", enables the trimming of reverse complement sequences of the provided primers. Default is `"true"`.
<!-- endif -->

    !!! techdetails "BBDuk Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_bbduk.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_bbduk.wdl) |
        | Software Source Code | [BBMap on SourceForge](https://sourceforge.net/projects/bbmap/) |
        | Software Documentation | [BBDuk Guide (archived)](https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) |
