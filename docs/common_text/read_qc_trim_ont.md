??? task "`read_QC_trim_ont`: Read Quality Trimming, Quantification, and Identification"

    `read_QC_trim_ont` is a sub-workflow that filters low-quality reads and trims low-quality regions of reads. It uses several tasks, described below.

<!-- if: theiacov|freyja -->
{{ include_md("common_text/ncbi_scrub_task.md", indent=4)}}

    ??? toggle "Read quality filtering"
      
        Read filtering is performed using `artic guppyplex` which performs a quality check by filtering the reads by length to remove chimeric reads.
<!-- endif -->

{{ include_md("common_text/kraken2_task.md", indent=4, replacements={"??? task": "??? toggle"}) }}
  
<!-- if: theiaprok -->
    !!! dna "A note on estimated genome length"

        By default, an estimated genome length is set to 5 Mb, which is around 0.7 Mb higher than the average bacterial genome length, according to the information collated [here](https://github.com/CDCgov/phoenix/blob/717d19c19338373fc0f89eba30757fe5cfb3e18a/assets/databases/NCBI_Assembly_stats_20240124.txt). This estimate can be overwritten by the user, and is used by `RASUSA`.
<!-- endif -->

    ??? toggle "`nanoplot`: Plotting and quantifying long-read sequencing data"

        Nanoplot is used for the determination of mean quality scores, read lengths, and number of reads. This task is run once with raw reads as input and once with clean reads as input. If QC has been performed correctly, you should expect **fewer** clean reads than raw reads.

<!-- if: theiaprok -->
    ??? toggle "Read subsampling"
        Samples are automatically randomly subsampled to 150X coverage using `RASUSA`.

    ??? toggle "Plasmid prediction"
        Plasmids are identified using replicon sequences used for typing from [PlasmidFinder](https://cge.food.dtu.dk/services/PlasmidFinder/).

    ??? toggle "Read filtering"
        Reads are filtered by length and quality using `nanoq`. By default, sequences with less than 500 basepairs and quality score lower than 10 are filtered out to improve assembly accuracy.
<!-- endif -->

    !!! techdetails "read_QC_trim_ont Technical Details"
        |  | Links |
        | --- | --- |
        | Sub-workflow | [wf_read_QC_trim_ont.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_read_QC_trim_ont.wdl) |
        | Tasks | [task_fastq_scan.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_fastq_scan.wdl)<br>[task_nanoplot.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_nanoplot.wdl)<br>[task_rasusa.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/task_rasusa.wdl)<br>[task_nanoq.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_nanoq.wdl)
        | Software Source Code | [fastq-scan](https://github.com/rpetit3/fastq-scan)<br>[NanoPlot](https://github.com/wdecoster/NanoPlot)<br>[RASUSA](https://github.com/mbhall88/rasusa)<br>[nanoq](https://github.com/esteinig/nanoq) |
        | Software Documentation | [NCBI Scrub on GitHub](<https://github.com/ncbi/sra-human-scrubber/blob/master/README.md>)<br>[NanoPlot documentation on GitHub](https://github.com/wdecoster/NanoPlot/tree/master?tab=readme-ov-file#readme)<br>[Rasusa on GitHub](https://github.com/mbhall88/rasusa)<br>[Nanoq on GitHub](https://github.com/esteinig/nanoq)<br>[Artic Pipeline ReadTheDocs](https://artic.readthedocs.io/en/latest/?badge=latest)<br>[Kraken2 on GitHub](https://github.com/DerrickWood/kraken2/wiki)
        | Original Publication(s) | [NanoPack2: population-scale evaluation of long-read sequencing data](https://academic.oup.com/bioinformatics/article/39/5/btad311/7160911)<br>[Rasusa: Randomly subsample sequencing reads to a specified coverage](https://doi.org/10.21105/joss.03941)<br>[Nanoq: ultra-fast quality control for nanopore reads](https://doi.org/10.21105/joss.02991)<br>[Improved metagenomic analysis with Kraken 2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0) |
