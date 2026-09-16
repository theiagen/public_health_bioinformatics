---
title: Task Fragment `ksnp`
fragment: true
---
<!-- if: ksnp3 -->
??? task "`kSNP3`: Phylogenetic Construction"
<!-- endif -->
<!-- if: ksnp4 -->
??? task "`kSNP4`: Phylogenetic Construction"
    !!! tip "Ensuring Disk Space for Large Datasets"
        For kSNP4, disk space requirements scale with the number and size of input samples. By default, 500 GB of disk space is available for the ksnp4 task. For a genome size of 12 Mb, we recommend budgeting ~1 GB per genome, ensuring that the final disk size is larger than the set size you’ve selected.

<!-- endif -->
    This workflow is run on a set of assembly files to produce both pan-genome and core-genome phylogenies. This also results in alignment files which are used by downstream tasks.

<!-- if: ksnp3 -->
    !!! techdetails "kSNP3 Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_ksnp3.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/task_ksnp3.wdl) |
        | Software Source Code | [kSNP on SourceForge via the WayBack Machine](https://web.archive.org/web/20221006234124/https://sourceforge.net/projects/ksnp/files/)  |
        | Software Documentation | [kSNP on SourceForge via the WayBack Machine](https://web.archive.org/web/20221006234124/https://sourceforge.net/projects/ksnp/files/) |
        | Original Publication(s) | [kSNP3.0: SNP detection and phylogenetic analysis of genomes without genome alignment or reference genome](https://doi.org/10.1093/bioinformatics/btv271) |

<!-- if: ksnp4 -->  
    !!! techdetails "kSNP4 Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_ksnp4.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/task_ksnp4.wdl) |
        | Software Source Code | [kSNP4 on SourceForge](https://sourceforge.net/projects/ksnp/files/)  |
        | Software Documentation | [kSNP4 on SourceForge](https://sourceforge.net/projects/ksnp/files/)  |
        | Original Publication(s) | [Building Phylogenetic Trees from Genome Sequences With kSNP4](https://doi.org/10.1093/molbev/msad235) |
<!-- endif -->
