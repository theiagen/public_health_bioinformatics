??? task "`abricate_abaum`: Plasmid Identification"
    *Acinetobacter* plasmids are not included in the PlasmidFinder database (see [the above section on Plasmid Identification](#plasmid-identification)). Instead, the [AcinetobacterPlasmidTyping](https://github.com/MehradHamidian/AcinetobacterPlasmidTyping) database contains variants of the plasmid *rep* gene for *A. baumannii* plasmid identification. When matched with >/= 95 % identity, this represents a typing scheme for *Acinetobacter baumannii* plasmids.

    The bioinformatics software for querying sample assemblies against the AcinetobacterPlasmidTyping database is [ABRicate](https://github.com/tseemann/abricate). By default, a 95% minimum identity threshold is set in order for successful classification.

    !!! techdetails "abricate_abaum Technical Details"

        |  | Links |
        | --- | --- |
        | Task | [task_abricate.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/drug_resistance/task_abricate.wdl) |
        | Software Source Code | [AcinetobacterPlasmidTyping Database on GitHub](https://github.com/MehradHamidian/AcinetobacterPlasmidTyping/tree/v1.0.0)<br>[ABRicate on GitHub](https://github.com/tseemann/abricate) |
        | Software Documentation | [AcinetobacterPlasmidTyping Database on GitHub](https://github.com/MehradHamidian/AcinetobacterPlasmidTyping/tree/v1.0.0)<br>[ABRicate on GitHub](https://github.com/tseemann/abricate) |
        | Original Publication(s) | _AcinetobacterPlasmidTyping database_: [Detection and Typing of Plasmids in *Acinetobacter baumannii* Using *rep* Genes Encoding Replication Initiation Proteins](https://journals.asm.org/doi/10.1128/spectrum.02478-22?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed) |
