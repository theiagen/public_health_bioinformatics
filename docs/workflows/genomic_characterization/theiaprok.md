# TheiaProk Workflow Series

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Name": "[**TheiaProk Workflow Series**](../workflows/genomic_characterization/theiaprok.md)"}, columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level", "Dockstore"]) }}

## TheiaProk Workflows

**The TheiaProk workflows are for the assembly, quality assessment, and characterization of bacterial genomes.** There are currently four TheiaProk workflows designed to accommodate different kinds of input data:

1. Illumina paired-end sequencing **(TheiaProk_Illumina_PE**)
2. Illumina single-end sequencing (**TheiaProk_Illumina_SE)**
3. ONT sequencing (**TheiaProk_ONT**)
4. Genome assemblies (**TheiaProk_FASTA**)

=== "TheiaProk_Illumina_PE"

    !!! caption "TheiaProk Illumina PE Workflow Diagram"
        ![TheiaProk Illumina PE Workflow Diagram](../../assets/figures/TheiaProk_v3.0.0.png)

=== "TheiaProk_ONT"

    !!! caption "TheiaProk ONT Workflow Diagram"
        ![TheiaProk ONT Workflow Diagram](../../assets/figures/TheiaProk_ONT.png)

All input reads are processed through "[core tasks](#core-tasks)" in the TheiaProk Illumina and ONT workflows. These undertake read trimming and assembly appropriate to the input data type. TheiaProk workflows subsequently launch default genome characterization modules for quality assessment, species identification, antimicrobial resistance gene detection, sequence typing, and more. **For some taxa identified, "[taxa-specific sub-workflows](#taxa-specific-tasks)" will be automatically activated, undertaking additional taxa-specific characterization steps**. When setting up each workflow, users may choose to use "optional tasks" as additions or alternatives to tasks run in the workflow by default.

### Inputs

!!! dna "Input Data"
    === "TheiaProk_Illumina_PE"

        The TheiaProk_Illumina_PE workflow takes in Illumina paired-end read data. Read file names should end with `.fastq` or `.fq`, with the optional addition of `.gz`. When possible, Theiagen recommends zipping files with [gzip](https://www.gnu.org/software/gzip/) before Terra uploads to minimize data upload time.

        By default, the workflow anticipates **2 x 150bp** reads (i.e. the input reads were generated using a 300-cycle sequencing kit). Modifications to the optional parameter for `trim_minlen` may be required to accommodate shorter read data, such as the 2 x 75bp reads generated using a 150-cycle sequencing kit.

    === "TheiaProk_Illumina_SE"

        TheiaProk_Illumina_SE takes in Illumina single-end reads. Read file names should end with `.fastq` or `.fq`, with the optional addition of `.gz`. Theiagen highly recommends zipping files with [gzip](https://www.gnu.org/software/gzip/) before uploading to Terra to minimize data upload time & save on storage costs.

        By default, the workflow anticipates **1 x 35 bp** reads  (i.e. the input reads were generated using a 70-cycle sequencing kit). Modifications to the optional parameter for `trim_minlen` may be required to accommodate longer read data.

    === "TheiaProk_ONT"

        The TheiaProk_ONT workflow takes in base-called ONT read data. Read file names should end with `.fastq` or `.fq`, with the optional addition of `.gz`. When possible, Theiagen recommends zipping files with [gzip](https://www.gnu.org/software/gzip/) before uploading to Terra to minimize data upload time.

        **The ONT sequencing kit and base-calling approach can produce substantial variability in the amount and quality of read data. Genome assemblies produced by the TheiaProk_ONT workflow must be quality assessed before reporting results.**

    === "TheiaProk_FASTA"

        The TheiaProk_FASTA workflow takes in assembly files in FASTA format.

!!! tip "Skip Characterization"
    Ever wanted to skip characterization? Now you can! Set the optional input `perform_characterization` to **`false`** to only generate an assembly and run assembly QC.

!!! caption ""
    === "TheiaProk_Illumina_PE"
        /// html | div[class="searchable-table"]

        {{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "TheiaProk_Illumina_PE"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"], indent=8) }}
        ///

    === "TheiaProk_Illumina_SE"
        /// html | div[class="searchable-table"]

        {{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "TheiaProk_Illumina_SE"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"], indent=8) }}
        ///

    === "TheiaProk_ONT"
        /// html | div[class="searchable-table"]

        {{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "TheiaProk_ONT"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"], indent=8) }}
        ///

    === "TheiaProk_FASTA"
        /// html | div[class="searchable-table"]

        {{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "TheiaProk_FASTA"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"], indent=8) }}
        ///

### Core Tasks

!!! dna ""
    These tasks are performed regardless of organism. They include tasks that are performed regardless of and specific for the input data type. They perform read trimming and assembly appropriate to the input data type.

{{ include_md("common_text/versioning_task.md") }}

#### Assembly Tasks

!!! caption ""
    === "TheiaProk_Illumina_PE"

{{ include_md("common_text/concatenate_illumina_lanes_task.md", indent=8) }}
{{ include_md("common_text/read_screen_task.md", condition="theiaprok", indent=8) }}

        !!! dna ""
            These tasks clean then assemble the reads into a _de novo_ assembly and assess the quality of the assembly.

{{ include_md("common_text/read_qc_trim_illumina_wf.md", condition="theiaprok", indent=8) }}
{{ include_md("common_text/cg_pipeline_task.md", indent=8) }}
{{ include_md("common_text/digger_denovo_wf.md", indent=8) }}

    === "TheiaProk_Illumina_SE"

{{ include_md("common_text/concatenate_illumina_lanes_task.md", indent=8) }}
{{ include_md("common_text/read_screen_task.md", condition="theiaprok", indent=8) }}

        !!! dna ""
            These tasks clean then assemble the reads into a _de novo_ assembly and assess the quality of the assembly.

{{ include_md("common_text/read_qc_trim_illumina_wf.md", condition="theiaprok", indent=8) }}
{{ include_md("common_text/cg_pipeline_task.md", indent=8) }}
{{ include_md("common_text/digger_denovo_wf.md", indent=8) }}

    === "TheiaProk_ONT"

{{ include_md("common_text/read_screen_task.md", condition="theiaprok", indent=8) }}

        !!! dna ""
            These tasks clean then assemble the reads into a _de novo_ assembly and assess the quality of the assembly.

{{ include_md("common_text/read_qc_trim_ont_wf.md", condition="theiaprok", indent=8) }}
{{ include_md("common_text/flye_denovo_wf.md", condition="theiaprok", indent=8) }}

    === "TheiaProk_FASTA"
        !!! dna ""
            Since this workflow requires FASTA files as input, no assembly or read trimming is performed, and the workflow proceeds directly to the "post-assembly tasks" section below.

#### Post-Assembly Tasks

!!! dna ""
    The following tasks are performed for all taxa, regardless of the input data type. They include quality control, assembly characterization, and taxonomic identification.

##### Quality Control

{{ include_md("common_text/quast_task.md") }}

{{ include_md("common_text/busco_task.md") }}

{{ include_md("common_text/qc_check_task.md", condition="theiaprok")}}

{{ include_md("common_text/arln_stats_task.md") }}

##### Taxonomic Assignment

{{ include_md("common_text/gambit_task.md") }}

{{ include_md("common_text/kmerfinder_task.md") }}

{{ include_md("common_text/mummer_ani_task.md") }}

##### AMR Characterization

{{ include_md("common_text/amrfinderplus_task.md") }}

{{ include_md("common_text/resfinder_task.md") }}

{{ include_md("common_text/gamma_task.md") }}

{{ include_md("common_text/abricate_bacterial_task.md") }}

{{ include_md("common_text/amr_search_task.md", condition="theiaprok") }}

##### Sequence Type

{{ include_md("common_text/ts_mlst_task.md") }}

##### Assembly Annotation

!!! tip ""
    By default, `genome_annotation` is set to `"prokka"`. To use `bakta` instead, set `genome_annotation` to `"bakta"`. These tools are mutually exclusive.

{{ include_md("common_text/prokka_task.md") }}

{{ include_md("common_text/bakta_task.md") }}

##### Plasmid Identification

{{ include_md("common_text/plasmidfinder_task.md") }}

### Taxa-Specific Tasks

!!! dna ""
    The TheiaProk workflows automatically activate taxa-specific sub-workflows after the identification of relevant taxa using `GAMBIT`. Alternatively, the user can provide the expected taxa in the `expected_taxon` workflow input to override the taxonomic assignment made by GAMBIT. Modules are launched for all TheiaProk workflows unless otherwise indicated.

    Please note that some modules require specific input data that may render it incompatible with every workflow. For example, ShigaTyper (a _Shigella_/EIEC serotyping tool) is not available for TheiaProk_FASTA as it requires read data as input. We have made a note next to each module to indicate which workflows are compatible with the module if there are restrictions.

{{ include_md("common_text/taxon_table_task.md") }}

??? toggle "_Acinetobacter baumannii_"
    ##### _Acinetobacter baumannii_ {% raw %} {#acinetobacter-baumannii} {% endraw %}

{{ include_md("common_text/kaptive_task.md", indent=4) }}

{{ include_md("common_text/abricate_abaum_task.md", indent=4) }}

??? toggle "_Escherichia_ or _Shigella_ spp."
    ##### _Escherichia_ or _Shigella_ spp. {% raw %} {#escherichia-or-shigella} {% endraw %}

    The *Escherichia* and *Shigella* genera are [difficult to differentiate as they do not comply with genomic definitions of genera and species](https://www.sciencedirect.com/science/article/abs/pii/S1286457902016374). Consequently, when either _Escherichia_ or _Shigella_ are identified by GAMBIT, all tools intended for these taxa are used. 

    !!! tip ""
        `SerotypeFinder` and `ECTyper` are intended for analysis of *E. coli*. Both tools are used as there are occasional discrepancies between the serotypes predicted. This primarily arises due to differences in the databases used by each tool.

{{ include_md("common_text/serotypefinder_task.md", indent=4)}}

{{ include_md("common_text/ectyper_task.md", indent=4)}}

    !!! tip ""
        `VirulenceFinder` identifies virulence genes in total or partial sequenced isolates of bacteria. Currently, only *E. coli* is supported in TheiaProk workflows. 

{{ include_md("common_text/virulencefinder_task.md", indent=4) }}

    !!! tip ""
        `ShigaTyper` and `ShigEiFinder` are intended for differentiation and serotype prediction for any *Shigella* species and Enteroinvasive *Escherichia coli* (EIEC). You can read about differences between these [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC517479/) and [here](https://www.microbiologyresearch.org/content/journal/micro/10.1099/00221287-144-9-2667). ShigEiFinder can be run using either the assembly (default) or reads. These tasks will report if the samples are neither *Shigella* nor EIEC.

{{ include_md("common_text/shigatyper_task.md", indent=4) }}

{{ include_md("common_text/shigeifinder_task.md", indent=4, condition="assembly") }}

{{ include_md("common_text/shigeifinder_task.md", indent=4, condition="reads") }}

    !!! tip ""
        `SonneiTyper` is run only when GAMBIT predicts the *S. sonnei* species. This is the most common *Shigella* species in the United States.

{{ include_md("common_text/sonneityper_task.md", indent=4)}}

    !!! tip ""
        **Shigella XDR prediction.** Please see [the documentation section above for ResFinder](#amr-characterization) for details regarding this taxa-specific analysis. 

{{ include_md("common_text/stxtyper_task.md", indent=4) }}

??? toggle "_Haemophilus influenzae_"
    ##### _Haemophilus influenzae_ {% raw %} {#haemophilus-influenzae} {% endraw %}

{{ include_md("common_text/hicap_task.md", indent=4) }}

??? toggle "_Klebsiella_ spp."
    ##### _Klebsiella_ spp. {% raw %} {#klebsiella} {% endraw %}

{{ include_md("common_text/kleborate_task.md", indent=4) }}

??? toggle "_Legionella pneumophila_"
    ##### _Legionella pneumophila_ {% raw %} {#legionella-pneumophila} {% endraw %}

{{ include_md("common_text/legsta_task.md", indent=4) }}

??? toggle "_Listeria monocytogenes_"
    ##### _Listeria monocytogenes_ {% raw %} {#listeria-monocytogenes} {% endraw %}

{{ include_md("common_text/lissero_task.md", indent=4) }}

??? toggle "_Mycobacterium tuberculosis_"
    ##### _Mycobacterium tuberculosis_ {% raw %} {#mycobacterium-tuberculosis} {% endraw %}

{{ include_md("common_text/clockwork_task.md", indent=4) }}

{{ include_md("common_text/tbprofiler_task.md", indent=4) }}

{{ include_md("common_text/tbp_parser_task.md", indent=4) }}

??? toggle "_Neisseria_ spp."
    ##### _Neisseria_ spp. {% raw %} {#neisseria} {% endraw %}

{{ include_md("common_text/ngmaster_task.md", indent=4) }}

{{ include_md("common_text/meningotype_task.md", indent=4) }}

??? toggle "_Pseudomonas aeruginosa_"
    ##### _Pseudomonas aeruginosa_ {% raw %} {#pseudomonas-aeruginosa} {% endraw %}

{{ include_md("common_text/pasty_task.md", indent=4) }}

??? toggle "_Salmonella_ spp."
    ##### _Salmonella_ spp. {% raw %} {#salmonella} {% endraw %}

    !!! tip ""
        Both SISTR and SeqSero2 are used for serotyping all *Salmonella* spp. Occasionally, the predicted serotypes may differ between SISTR and SeqSero2. When this occurs, differences are typically small and analogous, and are likely as a result of differing source databases. More information about Salmonella serovar nomenclature can be found [here](https://www.happykhan.com/posts/binfie-guide-serovar/).

{{ include_md("common_text/sistr_task.md", indent=4) }}

{{ include_md("common_text/seqsero2_task.md", indent=4) }}

{{ include_md("common_text/genotyphi_task.md", indent=4) }}

??? toggle "_Staphyloccocus aureus_"
    ##### _Staphyloccocus aureus_ {% raw %} {#staphyloccocus-aureus} {% endraw %}

{{ include_md("common_text/spatyper_task.md", indent=4) }}

{{ include_md("common_text/staphopia_sccmec_task.md", indent=4) }}

{{ include_md("common_text/agrvate_task.md", indent=4) }}

??? toggle "_Streptococcus pneumoniae_"
    ##### _Streptococcus pneumoniae_ {% raw %} {#streptococcus-pneumoniae} {% endraw %}

{{ include_md("common_text/poppunk_task.md", indent=4) }}

{{ include_md("common_text/seroba_task.md", indent=4) }}

{{ include_md("common_text/pbptyper_task.md", indent=4) }}

??? toggle "_Streptococcus pyogenes_"
    ##### _Streptococcus pyogenes_ {% raw %} {#streptococcus-pyogenes} {% endraw %}

{{ include_md("common_text/emmtypingtool_task.md", indent=4) }}

{{ include_md("common_text/emmtyper_task.md", indent=4) }}

??? toggle "_Vibrio_ spp."
    ##### _Vibrio_ spp. {% raw %} {#vibrio} {% endraw %}

{{ include_md("common_text/abricate_vibrio_task.md", indent=4) }}

{{ include_md("common_text/srst2_task.md", indent=4) }}

{{ include_md("common_text/vibecheck_task.md", indent=4) }}

### Outputs

=== "TheiaProk_Illumina_PE"
    /// html | div[class="searchable-table"]

    {{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "TheiaProk_Illumina_PE"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"], indent=4) }}

    ///

=== "TheiaProk_Illumina_SE"
    /// html | div[class="searchable-table"]

    {{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "TheiaProk_Illumina_SE"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"], indent=4) }}
    
    ///

=== "TheiaProk_ONT"
    /// html | div[class="searchable-table"]

    {{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "TheiaProk_FASTA"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"], indent=4) }}
    
    ///

=== "TheiaProk_FASTA"
    /// html | div[class="searchable-table"]

    {{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "TheiaProk_FASTA"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"], indent=4) }}

    ///
