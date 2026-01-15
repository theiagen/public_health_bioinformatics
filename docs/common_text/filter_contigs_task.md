<!-- if: flye -->
??? task "`Filter Contigs`: Filter contigs below a threshold length and remove homopolymer contigs"
    This task filters the created contigs based on a user-defined minimum length threshold (default of 1000) and eliminates homopolymer contigs (contigs of any length that consist of a single nucleotide).
<!-- endif -->
<!-- if: digger -->
??? task "`Filter Contigs`: Contig Quality Control"
    To **_de_**activate contig filtering, set `run_filter_contigs` to `false`.

    This task filters the created contigs based on a default minimum length threshold of 200 bp and a minimum coverage of 2.0. It also eliminates homopolymer contigs (contigs of any length that consist of a single nucleotide).

    Options are available to skip any of these filters by setting the respective parameters to `false`: `filter_contigs_skip_length_filter`, `filter_contigs_skip_coverage_filter`, and `filter_contigs_skip_homopolymer_filter`. The minimum length and coverage thresholds can be adjusted using the `filter_contigs_min_length` and `filter_contigs_min_coverage` parameters, respectively.
<!-- endif -->

    This ensures high-quality assemblies by retaining only contigs that meet specified criteria. Detailed metrics on contig counts and sequence lengths before and after filtering are provided in the output.

    !!! techdetails "Filter Contigs Technical Details" 
        |  | Links |
        | --- | --- |
        | WDL Task | [task_filter_contigs.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_filter_contigs.wdl) |
        