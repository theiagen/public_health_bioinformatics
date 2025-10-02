??? task "`ChroQueTas`: Antimicrobial Resistance Profiling"
    This task performs *in silico* antimicrobial resistance (AMR) mutation detection using [ChroQueTas](https://github.com/nmquijada/ChroQueTas) in conjunction with the [FungAMR](https://card.mcmaster.ca/fungamrhome) database.

    ChroQueTas extracts the putative expressed protein where point mutations are known to cause AMR, discards poorly aligned sequences via sequence similarity to known AMR genes via BLASTp, attempts to account for intron splicing and structural variation, and finally uses MAFFT to align and evaluate amino acid changes that may correspond to resistance relative to FungAMR genes.

    ChroQueTas is automatically called if a **compatible GAMBIT-predicted taxon or provided `expected_taxon` is identified**. Currently, 58 different species are supported.

    **Outputs:**

    - **AMR Stats**: Depicts the number of FungAMR and non-FungAMR mutations on a gene-by-gene (protein sequence) basis.
    - **AMR Summary**: Depicts the mutation resistance and support for potential drug resistance on a gene-by-gene (protein sequence) basis.
    - **Fungicide Resistance String**: Depicts detected AMR mutations AND the affected fungicides as a comma-delimited string in the following format: `"<PROTEIN>_<REFERENCE_AA><POSITION><MUTATION_AA>[<FUNGICIDE1><CONFIDENCE_SCORE>;...],..."`


    ??? toggle "Supported Species and Genes"
        The following table shows the species name and screened genes as of ChroQueTas v1.0.0. Species will be automatically determined based on the GAMBIT predicted taxon or provided `expected_taxon` name.

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
        
    ??? toggle "Confidence score for antimicrobial resistance"
        For each mutation there can be one or several antifungals for which the resistance has been associated with resistance in the past. For each antigungal, e.g. `Clotrimazole(3/NA)`, it appears two numbers (or an `NA`) between brackets separated by a `/`. These values represent the **confidence score**, that assesses the degree of evidence that supports their role in drug resistance. The number in the left side of the / is the **highest positive confidence score** while the one on the right side is the **highest negative confidence score**.

        - A **positive confidence score** denotes mutations reported to confer resistance
        - A **negative confidence score** relates to mutations reported in susceptible strains

        The following table shows the confidence scores and their rationale.

        | Confidence score    | Description                                                                                                                                                                                         |
        |---------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
        | 1                   | The mutation was created in the same susceptible species background and effect on resistance was confirmed by measurements.                                                                         |
        | 2                   | The mutation was created in the same gene but expressed in another species with potentially non endogenous level of expression and effect on resistance was confirmed by experimental measurements. |
        | 3                   | The effect of the mutation was quantified by bulk competition assays such as Deep Mutational Scanning.                                                                                              |
        | 4                   | The mutation was identified by experimental evolution where its high frequency in independent replicates suggests it causes resistance.                                                             |
        | 5                   | The deletion of the gene causes resistance.                                                                                                                                                         |
        | 6                   | The overexpression or duplication of the gene causes resistance.                                                                                                                                    |
        | 7                   | Significant association between the mutation and resistance in a population as determined by GWAS (or other population-based association such as QTL or classical genetics).                        |
        | 8                   | The mutation was identified in a 'natural' strain (e.g. a clinical isolate) that is resistant but without any further validations.                                                                  |
        | NA                  | Absence of a positive confidence score report                                                                                                                                                       |
        | Negative scores (-) | Any confidence score where the mutation was found not to cause resistance using the approach of the corresponding positive score.                                                                   |

        More information can be found in [ChroQueTas' documentation](https://github.com/nmquijada/ChroQueTas/wiki/Confidence-score-for-antimicrobial-resistance).

    !!! techdetails "ChroQueTas Technical Details"    
        |  | Links |
        | --- | --- |
        | Task | [task_chroquetas.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/drug_resistance/task_chroquetas.wdl) |
        | Software Source Code | [ChroQueTas on GitHub](https://github.com/nmquijada/ChroQueTas) |
        | Software Documentation | [ChroQueTas on GitHub](https://github.com/nmquijada/ChroQueTas) |
        | Original Publication(s) | [FungAMR: a comprehensive database for investigating fungal mutations associated with antimicrobial resistance](https://doi.org/10.1038/s41564-025-02084-7) |
