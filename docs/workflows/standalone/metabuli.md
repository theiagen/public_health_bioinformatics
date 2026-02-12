# Metabuli

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Name": "[**Metabuli**](../workflows/standalone/metabuli.md)"}, columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level", "Dockstore"]) }}

## Metabuli Workflow

**The Metabuli workflow assesses the taxonomic profile of raw sequencing data (FASTQ files).**

Metabuli is suitable for classifying short reads *AND* long reads, and optionally enables extracting reads from a specific NCBI taxon ID of interest. Metabuli uses a novel k-mer structure, called a "metamer", which incorporates both the DNA sequence for high specificity and amino acid conservation for sensitive homology detection.

The Metabuli_PHB workflow additionally includes read trimming software, Fastp (Illumina) and Porechop (ONT), for adapter trimming (recommended) and basic read preprocessing. 

!!! caption "Metabuli Workflow Diagram"
    ![Metabuli Workflow Diagram](../../assets/figures/Metabuli.png)

### Databases

!!! info  "Database selection"
    The Metabuli software is database-dependent and **taxonomic assignments are highly sensitive to the database used**. An appropriate database should contain the expected organism(s) (e.g. *Escherichia coli*) and other taxa that may be present in the reads (e.g. *Citrobacter freundii*, a common contaminant).

#### Suggested databases

<div class="searchable-table" markdown="1">

| Database name | Database Description | Suggested Applications | GCP URI (for usage in Terra) | Source | Database Size (GB) | Date of Last Update |
| --- | --- | --- | --- | --- | --- | --- |
| **viral** | RefSeq viral + human (T2T-CHM13v2.0) | Viral metagenomics | **`gs://theiagen-public-resources-rp/reference_data/databases/metabuli/refseq_virus-v223.tar.gz`** | <https://metabuli.steineggerlab.workers.dev/> | 4.0 | 2024/04/01 |
| **GTDB** | Prokaryote (Complete Genome/Chromosome, CheckM completeness > 90, and contamination <5) + human (T2T-CHM13v2.0) | Prokaryote metagenomics | **`gs://theiagen-public-resources-rp/reference_data/databases/metabuli/gtdb.tar.gz`** | <https://metabuli.steineggerlab.workers.dev/> | 68.8 | 2024/04/01 |

</div>

### Inputs

!!! caption ""
    === "Metabuli"
        /// html | div[class="searchable-table"]

        {{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "Metabuli"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"], indent=8) }}
        ///

??? dna "`taxon` input parameter"
    Inputting a `taxon` (NCBI taxon ID/name) will enable read extraction within the workflow. The input `taxon` will be standardized via querying the NCBI taxonomy hierarchy in the `ete4_identify` task. Additionally, a parent taxonomic `rank` (e.g. "genus", "family", "order", etc.) can be set in `ete4_identify` to extract reads at a higher taxonomic level relative to the input `taxon`.

??? dna "`illumina` input parameter"
    Setting `illumina` to "true" enables Illumina mode for single-end reads. Inputting a `read2` implicitly sets `illumina` to "true".

### Workflow Tasks

{{ include_md("common_text/ete4_identify_task.md") }}

{{ include_md("common_text/fastp_task.md", condition="metabuli") }}

{{ include_md("common_text/porechop_task.md", condition="metabuli") }}

{{ include_md("common_text/metabuli_task.md", condition="metabuli") }}

### Outputs

!!! caption ""
    === "Metabuli"
        /// html | div[class="searchable-table"]

        {{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "Metabuli"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"], indent=8) }}

        ///

#### Interpretation of results

The most important outputs of the Metabuli workflows are the `metabuli_report` files. These will include a breakdown of the number of sequences assigned to a particular taxon, and the percentage of reads assigned. [A complete description of the report format can be found here](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#standard-kraken-output-format).

When assessing the taxonomic identity of a single isolate's sequence, it is normal that a few reads are assigned to very closely rated taxa due to the shared sequence identity between them. "Very closely related taxa" may be genetically similar species in the same genus, or taxa with which the dominant species have undergone horizontal gene transfer. Unrelated taxa or a high abundance of these closely related taxa is indicative of contamination or sequencing of non-target taxa. Interpretation of the results is dependent on the biological context.

??? toggle "Example Metabuli report"
    Below is an example `metabuli_report` for a _Human immunodeficiency virus 1_ sample. Only the first 13 lines are included here since the rows near the bottom are <0.08% of the reads, which are likely human-derived contamination.

    From this report, we can see that ~98.78% of the reads were assigned at the species level (`species` in the 4th column) to "_Human immunodeficiency virus 1_". ~1.15% of the reads were unclassified, and the remaining <0.08% of reads are annoated as _Homo sapiens_ (not depicted).
    
    ```
     #clade_proportion	clade_count	taxon_count	rank	taxID	name
     1.1457	3045	3045	no rank	0	unclassified
     98.8543	262722	1	no rank	1	root
     98.7850	262538	0	superkingdom	10239	  Viruses
     98.7843	262536	0	clade	2559587	    Riboviria
     98.7843	262536	0	kingdom	2732397	      Pararnavirae
     98.7843	262536	0	phylum	2732409	        Artverviricota
     98.7843	262536	0	class	2732514	          Revtraviricetes
     98.7843	262536	0	order	2169561	            Ortervirales
     98.7843	262536	0	family	11632	              Retroviridae
     98.7843	262536	0	subfamily	327045	                Orthoretrovirinae
     98.7843	262536	0	genus	11646	                  Lentivirus
     **98.7843	262536	262536	species	11676	                    Human immunodeficiency virus 1**
    ```

#### Krona visualisation of Metabuli report

[Krona](https://github.com/marbl/Krona) produces an interactive report that allows hierarchical data, such as the one from Metabuli, to be explored with zooming, multi-layered pie charts. These pie charts are intuitive and highly responsive.

??? toggle "Example Krona report"

    Below is an example of the `krona_html` for a bacterial sample. Taxonomic rank is organised from the centre of the pie chart to the edge, with each slice representing the relative abundance of a given taxa in the sample.
    
    ![Example Krona Report](../../assets/figures/example_krona_report.png)

!!! techdetails "Metabuli Technical Details"
    |  | Links |
    | --- | --- |
    | Software Source Code | [Metabuli on GitHub](https://github.com/steineggerlab/metabuli/)  |
    | Software Documentation | <https://github.com/steineggerlab/Metabuli/blob/master/README.md> |
    | Original Publication(s) | [Metabuli: sensitive and specific metagenomic classification via joint analysis of amino acid and DNA](https://doi.org/10.1038/s41592-024-02273-y) |
