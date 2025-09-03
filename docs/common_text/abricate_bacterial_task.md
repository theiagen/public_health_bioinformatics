??? task "`ABRicate`: AMR Genotyping (optional)"

    The `abricate` module, if enabled, will run ABRicate with the database defined in `abricate_db` to perform mass screening of contigs for antimicrobial resistance or virulence genes. It comes bundled with multiple databases: NCBI, CARD, ARG-ANNOT, Resfinder, MEGARES, EcOH, PlasmidFinder, Ecoli_VF and VFDB. It only detects acquired resistance genes, **NOT** point mutations.

    By default, the "vfdb" database is used. The [virulence factor database (VFDB)](https://www.mgc.ac.cn/VFs/) is a comprehensive resource for bacterial pathogen virulence factors.

    !!! techdetails "ABRicate Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_abricate.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/drug_resistance/task_abricate.wdl) |
        | Software Source Code | [VFDB Database](https://www.mgc.ac.cn/VFs/main.htm)<br>[ABRicate on GitHub](https://github.com/tseemann/abricate) |
        | Software Documentation | [VFDB Database](https://www.mgc.ac.cn/VFs/main.htm)<br>[ABRicate on GitHub](https://github.com/tseemann/abricate) |
        | Original Publication(s) | _VFDB Database_: [VFDB 2019: a comparative pathogenomic platform with an interactive web interface](https://doi.org/10.1093/nar/gky1080) |
