# TheiaEuk Workflow Series

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibliity** | **Workflow Level** |
|---|---|---|---|---|
| [Genomic Characterization](../../workflows_overview/workflows_type.md/#genomic-characterization) | [Mycotics](../../workflows_overview/workflows_kingdom.md/#mycotics) | PHB v3.0.1 | Yes | Sample-level |

## TheiaEuk Workflows

**The TheiaEuk_Illumina_PE workflow is for the assembly, quality assessment, and characterization of fungal genomes.** It is designed to accept Illumina paired-end sequencing data as the primary input. **It is currently intended only for ==haploid== fungal genomes like _Candidozyma auris_.** Analyzing diploid genomes using TheiaEuk should be attempted only with expert attention to the resulting genome quality.

All input reads are processed through "core tasks" in each workflow. The core tasks include raw read quality assessment, read cleaning (quality trimming and adapter removal), de novo assembly, assembly quality assessment, and species taxon identification. For some taxa identified, taxa-specific sub-workflows will be automatically activated, undertaking additional taxa-specific characterization steps, including clade-typing and/or antifungal resistance detection.

!!! caption "TheiaEuk Workflow Diagram"
    ![TheiaEuk Workflow Diagram](../../assets/figures/TheiaEuk_Illumina_PHB_2025123.png){width=75%}

### Inputs

!!! info "Input read data"

    The TheiaEuk_Illumina_PE workflow takes in Illumina paired-end read data. Read file names should end with `.fastq` or `.fq`, with the optional addition of `.gz`. When possible, Theiagen recommends zipping files with [gzip](https://www.gnu.org/software/gzip/) prior to Terra upload to minimize data upload time.

    By default, the workflow anticipates 2 x 150bp reads (i.e. the input reads were generated using a 300-cycle sequencing kit). Modifications to the optional parameter for `trim_minlen` may be required to accommodate shorter read data, such as the 2 x 75bp reads generated using a 150-cycle sequencing kit.

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="TheiaEuk_Illumina_PE", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Workflow Tasks

All input reads are processed through "core tasks" in the TheiaEuk workflows. These undertake read trimming and assembly appropriate to the input data type, currently only Illumina paired-end data. TheiaEuk workflow subsequently launch default genome characterization modules for quality assessment, and additional taxa-specific characterization steps. When setting up the workflow, users may choose to use "optional tasks" or alternatives to tasks run in the workflow by default.

#### Core tasks

!!! tip ""
    These tasks are performed regardless of organism. They perform read trimming and various quality control steps.

{{ include_md("common_text/versioning_task.md") }}
{{ include_md("common_text/read_screen_task.md", condition="theiaeuk") }}

??? task "`Rasusa`: Read subsampling (optional, on by default)"

    The Rasusa task performs subsampling of the raw reads. By default, this task will subsample reads to a depth of 150X using the estimated genome length produced during the preceding raw read screen. The user can prevent the task from being launched by setting the `call_rasusa`variable to false. 

    The user can also provide an estimated genome length for the task to use for subsampling using the `genome_size` variable. In addition, the read depth can be modified using the `subsample_coverage` variable.
        
    !!! techdetails "Rasusa Technical Details"

        |  | Links |
        | --- | --- |
        | Task | [task_rasusa.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/task_rasusa.wdl) |
        | Software Source Code | [Rasusa on GitHub](https://github.com/mbhall88/rasusa) |
        | Software Documentation | [Rasusa on GitHub](https://github.com/mbhall88/rasusa) |
        | Original Publication(s) | [Rasusa: Randomly subsample sequencing reads to a specified coverage](https://doi.org/10.21105/joss.03941) |

{{ include_md("common_text/read_qc_trim_illumina.md", condition="theiaprok") }}

#### Assembly tasks

!!! tip ""
    These tasks assemble the reads into a _de novo_ assembly and assess the quality of the assembly.

{{ include_md("common_text/shovill_task.md") }}

??? task "`QUAST`: Assembly Quality Assessment"

    [`QUAST`](https://github.com/ablab/quast) (**QU**ality **AS**sessment **T**ool) evaluates genome assemblies by computing several metrics that describe the assembly quality, including the total number of bases in the assembly, the length of the largest contig in the assembly, and the assembly percentage GC content.

    !!! techdetails "QUAST Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_quast.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_quast.wdl) |
        | Software Source Code | [QUAST on GitHub](https://github.com/ablab/quast) |
        | Software Documentation | https://quast.sourceforge.net/docs/manual.html |
        | Orginal publication | [QUAST: quality assessment tool for genome assemblies](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3624806/) |

{{ include_md("common_text/cg_pipeline_task.md") }}

#### Organism-agnostic characterization

!!! tip ""
    These tasks are performed regardless of the organism and provide quality control and taxonomic assignment.

{{ include_md("common_text/gambit_task.md") }}

??? task "`BUSCO`: Assembly Quality Assessment"

    BUSCO (**B**enchmarking **U**niversal **S**ingle-**C**opy **O**rthologue) attempts to quantify the completeness and contamination of an assembly to generate quality assessment metrics. It uses taxa-specific databases containing genes that are all expected to occur in the given taxa, each in a single copy. BUSCO examines the presence or absence of these genes, whether they are fragmented, and whether they are duplicated (suggestive that additional copies came from contaminants).

    **BUSCO notation** 
    
    Here is an example of BUSCO notation: `C:99.1%[S:98.9%,D:0.2%],F:0.0%,M:0.9%,n:440`. There are several abbreviations used in this output:
    
    - Complete (C) - genes are considered "complete" when their lengths are within two standard deviations of the BUSCO group mean length.
    - Single-copy (S) - genes that are complete and have only one copy.
    - Duplicated (D) - genes that are complete and have more than one copy.
    - Fragmented (F) - genes that are only partially recovered.
    - Missing (M) - genes that were not recovered at all.
    - Number of genes examined (n) - the number of genes examined.
    
    A high equity assembly will use the appropriate database for the taxa, have high complete (C) and single-copy (S) percentages, and low duplicated (D), fragmented (F) and missing (M) percentages. 
  
    !!! techdetails "BUSCO Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_busco.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/advanced_metrics/task_busco.wdl) |
        | Software Source Code | [BUSCO on GitLab](https://gitlab.com/ezlab/busco) |
        | Software Documentation | https://busco.ezlab.org/ |
        | Orginal publication | [BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs](https://academic.oup.com/bioinformatics/article/31/19/3210/211866) |

{{ include_md("common_text/qc_check_task.md", condition="theiaeuk")}}

#### Organism-specific characterization

!!! tip ""
    The TheiaEuk workflow automatically activates taxa-specific tasks after identification of the relevant taxa using `GAMBIT`. Many of these taxa-specific tasks do not require any additional inputs from the user.

??? toggle "_Candidozyma auris_ (also known as _Candida auris_)"
    Two tools are deployed when _Candidozyma auris_/_Candida auris_ is  identified.

    ??? task "Cladetyping: clade determination"
        GAMBIT is used to determine the clade of the specimen by comparing the sequence to five clade-specific reference files. The output of the clade typing task will be used to specify the reference genome for the antifungal resistance detection tool.

        ??? toggle "Default reference genomes used for clade typing and antimicrobial resistance gene detection of _C. auris_"
            | Clade | Genome Accession | Assembly Name | Strain | NCBI Submitter | Included mutations in AMR genes (not comprehensive) |
            | --- | --- | --- | --- | --- | --- |
            | _Candidozyma auris_ Clade I | GCA_002759435.2 | Cand_auris_B8441_V2 | B8441 | Centers for Disease Control and Prevention |  |
            | _Candidozyma auris_ Clade II | GCA_003013715.2 | ASM301371v2 | B11220 | Centers for Disease Control and Prevention |  |
            | _Candidozyma auris_ Clade III | GCA_002775015.1 | Cand_auris_B11221_V1 | B11221 | Centers for Disease Control and Prevention | _ERG11_ V125A/F126L |
            | _Candidozyma auris_ Clade IV | GCA_003014415.1 | Cand_auris_B11243 | B11243 | Centers for Disease Control and Prevention | _ERG11_ Y132F |
            | _Candidozyma auris_ Clade V | GCA_016809505.1 | ASM1680950v1 | IFRC2087 | Centers for Disease Control and Prevention |  |

        !!! techdetails "Cladetyping Technical Details"
            |  | Links |
            | --- | --- |
            | Task | [task_cauris_cladetyping.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/candida/task_cauris_cladetyper.wdl) |
            | Software Source Code | [GAMBIT on GitHub](https://github.com/jlumpe/gambit) |
            | Software Documentation | [GAMBIT Overview](https://theiagen.notion.site/GAMBIT-7c1376b861d0486abfbc316480046bdc?pvs=4)
            | Original Publication(s) | [GAMBIT (Genomic Approximation Method for Bacterial Identification and Tracking): A methodology to rapidly leverage whole genome sequencing of bacterial isolates for clinical identification](https://doi.org/10.1371/journal.pone.0277575)<br> [TheiaEuk: a species-agnostic bioinformatics workflow for fungal genomic characterization](https://doi.org/10.3389/fpubh.2023.1198213) |

    ??? task "Snippy Variants: antifungal resistance detection"
        To detect mutations that may confer antifungal resistance, `Snippy` is used to find all variants relative to the clade-specific reference, then these variants are queried for product names associated with resistance.
    
        The genes in which there are known resistance-conferring mutations for this pathogen are:

        - FKS1
        - ERG11 (lanosterol 14-alpha demethylase)
        - FUR1 (uracil phosphoribosyltransferase)

        We query `Snippy` results to see if any mutations were identified in those genes. By default, we automatically check for the following loci (which can be overwritten by the user). You will find the mutations next to the locus tag in the `theiaeuk_snippy_variants_hits` column corresponding gene name (see below):

        | **TheiaEuk Search Term** | **Corresponding Gene Name** |
        |---|---|
        | B9J08_005340 | ERG6 |
        | B9J08_000401 | FLO8 |
        | B9J08_005343 | Hypothetical protein (PSK74852) |
        | B9J08_003102 | MEC3 |
        | B9J08_003737 | ERG3 |
        | lanosterol.14-alpha.demethylase | ERG11 |
        | uracil.phosphoribosyltransferase | FUR1 |
        | FKS1 | FKS1 |    

        For example, one sample may have the following output for the `theiaeuk_snippy_variants_hits` column:

        ```plaintext
        lanosterol.14-alpha.demethylase: lanosterol 14-alpha demethylase (missense_variant c.428A>G p.Lys143Arg; C:266 T:0),B9J08_000401: hypothetical protein (stop_gained c.424C>T p.Gln142*; A:70 G:0)
        ```

        Based on this, we can tell that ERG11 has a missense variant at position 143 (Lysine to Arginine) and B9J08_000401 (which is FLO8) has a stop-gained variant at position 142 (Glutamine to Stop).

        ??? toggle "Known resistance-conferring mutations for _Candidozyma auris_"
            Mutations in these genes that are known to confer resistance are shown below

            | **Organism** | **Found in** | **Gene name** | **Gene locus** | **AA mutation** | **Drug** | **Reference** |
            | --- | --- | --- | --- | --- | --- | --- |
            | **Candidozyma auris** | **Human** | **ERG11** |  | **Y132F** | **Fluconazole** | [Simultaneous Emergence of Multidrug-Resistant _Candida auris_ on 3 Continents Confirmed by Whole-Genome Sequencing and Epidemiological Analyses](https://academic.oup.com/cid/article/64/2/134/2706620/Simultaneous-Emergence-of-Multidrug-Resistant) |
            | **Candidozyma auris** | **Human** | **ERG11** |  | **K143R** | **Fluconazole** | [Simultaneous Emergence of Multidrug-Resistant _Candida auris_ on 3 Continents Confirmed by Whole-Genome Sequencing and Epidemiological Analyses](https://academic.oup.com/cid/article/64/2/134/2706620/Simultaneous-Emergence-of-Multidrug-Resistant) |
            | **Candidozyma auris** | **Human** | **ERG11** |  | **F126T** | **Fluconazole** | [Simultaneous Emergence of Multidrug-Resistant _Candida auris_ on 3 Continents Confirmed by Whole-Genome Sequencing and Epidemiological Analyses](https://academic.oup.com/cid/article/64/2/134/2706620/Simultaneous-Emergence-of-Multidrug-Resistant) |
            | **Candidozyma auris** | **Human** | **FKS1** |  | **S639P** | **Micafungin**  | [Activity of CD101, a long-acting echinocandin, against clinical isolates of Candida auris](https://www.sciencedirect.com/science/article/pii/S0732889317303498) |
            | **Candidozyma auris** | **Human** | **FKS1** |  | **S639P** | **Caspofungin** | [Activity of CD101, a long-acting echinocandin, against clinical isolates of Candida auris](https://www.sciencedirect.com/science/article/pii/S0732889317303498) |
            | **Candidozyma auris** | **Human** | **FKS1** |  | **S639P** | **Anidulafungin** | [Activity of CD101, a long-acting echinocandin, against clinical isolates of Candida auris](https://www.sciencedirect.com/science/article/pii/S0732889317303498) |
            | **Candidozyma auris** | **Human** | **FKS1** |  | **S639F** | **Micafungin** | [A multicentre study of antifungal susceptibility patterns among 350 _Candida auris_ isolates (2009–17) in India: role of the ERG11 and FKS1 genes in azole and echinocandin resistance](https://academic.oup.com/jac/advance-article/doi/10.1093/jac/dkx480/4794718) |
            | **Candidozyma auris** | **Human** | **FKS1** |  | **S639F** | **Caspofungin** | [A multicentre study of antifungal susceptibility patterns among 350 _Candida auris_ isolates (2009–17) in India: role of the ERG11 and FKS1 genes in azole and echinocandin resistance](https://academic.oup.com/jac/advance-article/doi/10.1093/jac/dkx480/4794718) |
            | **Candidozyma auris** | **Human** | **FKS1** |  | **S639F** | **Anidulafungin** | [A multicentre study of antifungal susceptibility patterns among 350 _Candida auris_ isolates (2009–17) in India: role of the ERG11 and FKS1 genes in azole and echinocandin resistance](https://academic.oup.com/jac/advance-article/doi/10.1093/jac/dkx480/4794718) |
            | **Candidozyma auris** | **Human** | **FUR1** | **CAMJ_004922** | **F211I** | **5-flucytosine** | [Genomic epidemiology of the UK outbreak of the emerging human fungal pathogen Candida auris](https://doi.org/10.1038/s41426-018-0045-x) |

        !!! techdetails "Snippy Variants Technical Details"
            |  | Links |
            | --- | --- |
            | Task | [task_snippy_variants.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snippy_variants.wdl)<br>[task_snippy_gene_query.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snippy_gene_query.wdl) |
            | Software Source Code | [Snippy on GitHub](https://github.com/tseemann/snippy) |
            | Software Documentation | [Snippy on GitHub](https://github.com/tseemann/snippy) |
            
??? toggle "_Candida albicans_"
    When this species is detected by the taxon ID tool, an antifungal resistance detection task is deployed.

    ??? task "Snippy Variants: antifungal resistance detection"
        To detect mutations that may confer antifungal resistance, `Snippy` is used to find all variants relative to the clade-specific reference, and these variants are queried for product names associated with resistance.

        The genes in which there are known resistance-conferring mutations for this pathogen are:

        - ERG11
        - GCS1 (FKS1)
        - FUR1
        - RTA2

        We query `Snippy` results to see if any mutations were identified in those genes. By default, we automatically check for the following loci (which can be overwritten by the user). You will find the mutations next to the locus tag in the `theiaeuk_snippy_variants_hits` column corresponding gene name (see below):

        | **TheiaEuk Search Term** | **Corresponding Gene Name** |
        |---|---|
        | ERG11 | ERG11 |
        | GCS1 | FKS1 |
        | FUR1 | FUR1 |
        | RTA2 | RTA2 |

        !!! techdetails "Snippy Variants Technical Details"
            |  | Links |
            | --- | --- |
            | Task | [task_snippy_variants.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snippy_variants.wdl)<br>[task_snippy_gene_query.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snippy_gene_query.wdl) |
            | Software Source Code | [Snippy on GitHub](https://github.com/tseemann/snippy) |
            | Software Documentation | [Snippy on GitHub](https://github.com/tseemann/snippy) |

??? toggle "_Aspergillus fumigatus_"
    When this species is detected by the taxon ID tool an antifungal resistance detection task is deployed.

    ??? task "Snippy Variants: antifungal resistance detection"
        To detect mutations that may confer antifungal resistance, `Snippy` is used to find all variants relative to the clade-specific reference, and these variants are queried for product names associated with resistance.

        The genes in which there are known resistance-conferring mutations for this pathogen are:

        - Cyp51A
        - HapE
        - COX10 (AFUA_4G08340)
 
        We query `Snippy` results to see if any mutations were identified in those genes. By default, we automatically check for the following loci (which can be overwritten by the user). You will find the mutations next to the locus tag in the `theiaeuk_snippy_variants_hits` column corresponding gene name (see below):

        | **TheiaEuk Search Term** | **Corresponding Gene Name** |
        |---|---|
        | Cyp51A | Cyp51A |
        | HapE | HapE |
        | AFUA_4G08340 | COX10 |

        !!! techdetails "Snippy Variants Technical Details"
            |  | Links |
            | --- | --- |
            | Task | [task_snippy_variants.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snippy_variants.wdl)<br>[task_snippy_gene_query.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snippy_gene_query.wdl) |
            | Software Source Code | [Snippy on GitHub](https://github.com/tseemann/snippy) |
            | Software Documentation | [Snippy on GitHub](https://github.com/tseemann/snippy) |

??? toggle "_Cryptococcus neoformans_"
    When this species is detected by the taxon ID tool an antifungal resistance detection task is deployed.

    ??? task "Snippy Variants: antifungal resistance detection"
        To detect mutations that may confer antifungal resistance, `Snippy` is used to find all variants relative to the clade-specific reference, and these variants are queried for product names associated with resistance.

        The genes in which there are known resistance-conferring mutations for this pathogen are:

        - ERG11 (CNA00300)
        
        We query `Snippy` results to see if any mutations were identified in those genes. By default, we automatically check for the following loci (which can be overwritten by the user). You will find the mutations next to the locus tag in the `theiaeuk_snippy_variants_hits` column corresponding gene name (see below):

        | **TheiaEuk Search Term** | **Corresponding Gene Name** |
        |---|---|
        | CNA00300 | ERG11 |
    
        !!! techdetails "Snippy Variants Technical Details"
            |  | Links |
            | --- | --- |
            | Task | [task_snippy_variants.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snippy_variants.wdl)<br>[task_snippy_gene_query.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_snippy_gene_query.wdl) |
            | Software Source Code | [Snippy on GitHub](https://github.com/tseemann/snippy) |
            | Software Documentation | [Snippy on GitHub](https://github.com/tseemann/snippy) |

### Outputs

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/tables/all_outputs.tsv", input_table=False, filter_column="Workflow", filter_values="TheiaEuk_Illumina_PE", columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///
