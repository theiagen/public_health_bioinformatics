<!-- if: long_read_flags|only_map_ont -->
??? task "`minimap2`: Read Alignment Details"
<!-- endif -->

    `minimap2` is a popular aligner that is used to align reads (or assemblies) to an assembly file. In minimap2, "modes" are a group of preset options.

<!-- if: long_read_flags -->
    The mode used in this task is `map-ont` with additional long-read-specific parameters (the `-L --cs --MD` flags) to align ONT reads to the reference genome. These specialized parameters are essential for proper handling of long read error profiles, generation of detailed alignment information, and improved mapping accuracy for long reads.

    `map-ont` is the default mode for long reads and it indicates that long reads of ~10% error rates should be aligned to the reference genome. The output file is in SAM format.
<!-- endif -->

<!-- if: only_map_ont -->
    The mode used in this task is `map-ont` which is the default mode for long reads and indicates that long reads of ~10% error rates should be aligned to the reference genome. The output file is in SAM format.
<!-- endif -->

<!-- if: asm20_mode -->
    The mode used in this task is `asm20` which is intended for "long assembly to reference mapping". The `asm20` mode indicates the following parameters should be used: `-k19 -w10 -U50,500 --rmq -r100k -g10k -A1 -B4 -O6,26 -E2,1 -s200 -z200 -N50`. The output file is in PAF format.
<!-- endif -->

<!-- if: sr_mode -->
    The mode used in this task is `sr` which is intended for "short single-end reads without splicing". The `sr` mode indicates the following parameters should be used: `-k21 -w11 --sr --frag=yes -A2 -B8 -O12,32 -E2,1 -b0 -r100 -p.5 -N20 -f1000,5000 -n2 -m20 -s40 -g100 -2K50m --heap-sort=yes --secondary=no`. The output file is in SAM format.
<!-- endif -->

    For more information regarding modes and the available options for `minimap2`, please see the [minimap2 manpage](https://lh3.github.io/minimap2/minimap2.html)

    !!! techdetails "minimap2 Technical Details"
        | | Links |
        |---|---|
        | Task | [task_minimap2.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/alignment/task_minimap2.wdl) |
        | Software Source Code | [minimap2 on GitHub](https://github.com/lh3/minimap2) |
        | Software Documentation | [minimap2](https://lh3.github.io/minimap2) |
        | Original Publication(s) | [Minimap2: pairwise alignment for nucleotide sequences](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778) |
