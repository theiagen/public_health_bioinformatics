??? task "`chroquetas`: Antimicrobial Resistance Profiling"
    This task performs *in silico* antimicrobial resistance (AMR) mutation detection using ChroQueTas in conjunction with the [FungAMR](https://card.mcmaster.ca/fungamrhome) database.

    ChroQueTas extracts the putative expressed protein where point mutations are known to cause AMR, discards poorly aligned sequences via sequence similarity to known AMR genes via BLASTp, attempts to account for intron splicing and structural variation, and finally uses MAFFT to align and evaluate amino acid changes that may correspond to resistance relative to FungAMR genes.

    ChroQueTas is automatically called if a compatible GAMBIT-predicted taxon is identified.

    **Outputs:**

    - **AMR Stats**: Depicts the number of FungAMR and non-FungAMR mutations on a gene-by-gene (protein sequence) basis.
    - **AMR Summary**: Depicts the mutation resistance and support for potential drug resistance on a gene-by-gene (protein sequence) basis. To interpret confidence scores, please refer to the corresponding [ChroQueTas wiki](https://github.com/nmquijada/ChroQueTas/wiki/Confidence-score-for-antimicrobial-resistance) section.
    - **AMR String**: Depicts detected AMR mutations as a comma-delimited string in the following format: "<PROTEIN>_<REFERNCE_AA><POSITION><MUTATION_AA>,..."


    ??? toggle "Supported Species and Genes"
        The following table shows the species name and screened genes as of ChroQueTas v1.0.0. Species will be automatically determined based on the GAMBIT predicted taxon.

    | Species                          | Antimicrobial Resistance Genes                                      |
    |-----------------------------------|---------------------------------------------------------------------|
    | *Arthroderma vanbreuseghemii*     | Squalene epoxidase                                                 |
    | *Aspergillus awamori*             | Cyp51                                                              |
    | *Aspergillus fumigatus*           | Beta-tubulin, Cytochrome b, Fks, HapE, Hmg1, Hmg2, Sdh             |
    | *Aspergillus niger*               | Cyp51                                                              |
    | *Aspergillus terreus*             | Cyp51                                                              |
    | *Aspergillus tubingensis*         | Cyp51                                                              |
    | *Beauveria bassiana*              | Beta-tubulin                                                       |
    | *Blumeria graminis*               | Cyp51                                                              |
    | *Blumeria graminis hordei*        | Cyp51                                                              |
    | *Blumeria graminis tritici*       | Cytochrome b                                                       |
    | *Botrytis cinerea*                | Beta-tubulin, Bos1, Cytochrome b, OS1                              |
    | *Candida albicans*                | Cap1, Cdr1, Cyp51, Erg3, Fur1, Mrr1, Mrr2, Rta2, Squalene epoxidase, Tac1, Upc2 |
    | *Candida dubliniensis*            | Cyp51, Erg3, Fcy1, Fks                                             |
    | *Candida metapsilosis*            | Fks                                                                |
    | *Candida orthopsilosis*           | Cyp51, Fks                                                         |
    | *Candida parapsilosis*            | Cyp51, Erg3, Fks, Mrr1                                             |
    | *Candida tropicalis*              | Ben1, Cyp51, Erg3, Fks                                             |
    | *Candidozyma auris*               | Cyp51, Erg10, Erg12, Erg3, Erg6, Fks, Fur1, Hmg1, Ncp1, Tac1       |
    | *Cercospora beticola*             | Beta-tubulin, Cytochrome b                                         |
    | *Clavispora lusitaniae*           | Fcy1, Fcy2, Fks                                                    |
    | *Colletotrichum acutatum*         | Cyp51, Cytochrome b                                                |
    | *Colletotrichum fructicola*       | Beta-tubulin, Cytochrome b                                         |
    | *Colletotrichum gloeosporioides*  | Beta-tubulin                                                       |
    | *Colletotrichum graminicola*      | Cytochrome b                                                       |
    | *Colletotrichum siamense*         | Cytochrome b                                                       |
    | *Corynespora cassiicola*          | Cytochrome b                                                       |
    | *Cryptococcus deuterogattii*      | Uxs1                                                               |
    | *Cryptococcus neoformans*         | Bck1, CNAG 04784, CNBG 2198, Cyp51, Dpb2, Fcy1, Fcy2, Fur1, Ran1, Ura6 |
    | *Erysiphe necator*                | Sdh                                                                |
    | *Fusarium fujikuroi*              | Cyp51                                                              |
    | *Fusarium graminearum*            | Beta-tubulin, Myosin-1, Myosin-5                                   |
    | *Fusarium incarnatum*             | Myosin-5                                                           |
    | *Fusarium pseudograminearum*      | Sdh                                                                |
    | *Histoplasma capsulatum*          | Cyp51                                                              |
    | *Kluyveromyces marxianus*         | Fks                                                                |
    | *Magnaporthe grisea*              | Cytochrome b, Sdh                                                  |
    | *Monilinia fructicola*            | Beta-tubulin                                                       |
    | *Mycosphaerella fijiensis*        | Cyp51, Cytochrome b                                                |
    | *Nakaseomyces glabratus*          | Ben1, Cdc55, Cdc6, Cdr1, Cyp51, Dot6, Erg2, Erg3, Erg4, Erg6, Erg7, Erg8, Fcy1, Fen1, Fur1, Gph1, Mohi, Mrpl11, Msh2, Pdr1, Pdr15, Qdr2, Sui2, Tcb |
    | *Passalora fulva*                 | Beta-tubulin                                                       |
    | *Phaeosphaeria nodorum*           | Cyp51                                                              |
    | *Phakopsora pachyrhizi*           | Cyp51, Cytochrome b                                                |
    | *Pichia kudriavzevii*             | Cyp51, Fks                                                         |
    | *Pneumocystis jirovecii*          | Dhfr                                                               |
    | *Puccinia triticina*              | Cyp51                                                              |
    | *Pyrenopeziza brassicae*          | Beta-tubulin                                                       |
    | *Saccharomyces cerevisiae*        | Act1, Aft1, Any1, Arg2, Atp1, Atp2, Bck1, Bcy1, Beta-tubulin, Bul1, Cbk1, Cdc43, Cdc60, Cdr1, Csf1, Csg2, Cyp51, Cyr1, Cytochrome b, Dcp2, Dhh1, Elo2, Erg12, Erg20, Erg25, Erg3, Erg6, Erg7, Erg9, Fcy1, Fcy21, Fpr1, Fur1, Gfa1, Hap1, Hem1, Hrd3, Hym1, Inp53, Lem3, Pde2, Pdr11, Pma1, Pms1, Pol2, Rav2, Rox1, Rpl28, Rsp5, Scd5, Sip3, Sit4, Sur1, Swh1, Tao3, Tfc1, Top1, Top2, Tsc11, Tup1, Ura6, Utp18, Vma1, Vma16, Vma9, Yrr1 |
    | *Saccharomyces paradoxus*         | Cdr1, Erg2, Erg6, Pdr10, Upc2, Yrr1                                |
    | *Tapesia acuformis*               | Cyp51                                                              |
    | *Tapesia yallundae*               | Cyp51                                                              |
    | *Trichophyton indotineae*         | Squalene epoxidase                                                 |
    | *Trichophyton interdigitale*      | Squalene epoxidase                                                 |
    | *Trichophyton mentagrophytes*     | Squalene epoxidase                                                 |
    | *Trichophyton rubrum*             | Squalene epoxidase                                                 |
    | *Ustilaginoidea virens*           | Cyp51                                                              |
    | *Venturia inaequalis*             | Beta-tubulin, Cytochrome b                                         |
    | *Zymoseptoria tritici*            | Cyp51, Cytochrome b                                                |
    

    !!! techdetails "ChroQueTas Technical Details"    
        |  | Links |
        | --- | --- |
        | Task | [task_chroquetas.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/drug_resistance/task_chroquetas.wdl) |
        | Software Source Code | [ChroQueTas on GitHub](https://github.com/nmquijada/ChroQueTas) |
        | Software Documentation | [ChroQueTas on GitHub](https://github.com/nmquijada/ChroQueTas) |
        | Original Publication(s) | [FungAMR: a comprehensive database for investigating fungal mutations associated with antimicrobial resistance](https://doi.org/10.1038/s41564-025-02084-7) |
