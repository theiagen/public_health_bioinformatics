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
| **viral** | RefSeq viral | Viral metagenomics | **`gs://theiagen-public-resources-rp/reference_data/databases/metabuli/refseq_virus-v223.tar.gz`** | <https://metabuli.steineggerlab.workers.dev/> | 4.0 | 2024/04/01 |
| **GTDB** | Prokaryote (Complete Genome/Chromosome, CheckM completeness > 90 and contamination <5) + human (T2T-CHM13v2.0) | **`gs://theiagen-public-resources-rp/reference_data/databases/metabuli/gtdb.tar.gz`** | <https://metabuli.steineggerlab.workers.dev/> | 68.8 | 2024/04/01 |

</div>

### Inputs

!!! caption ""
    === "Metabuli"
        /// html | div[class="searchable-table"]

        {{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "Metabuli"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"], indent=8) }}
        ///

??? dna "`illumina` input parameter
   Setting `illumina` to "true" enables Illumina mode for single-end reads. Inputting a `read2` implicitly sets `illumina` to "true".

### Workflow Tasks

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
    Below is an example `Metabuli_report` for a _Klebsiella pneumoniae_ sample. Only the first 30 lines are included here since rows near the bottom are often spurious results with only a few reads assigned to a non-target organism.

    From this report, we can see that 84.35 % of the reads were assigned at the species level (`S` in the 4th column) to "_Klebsiella pneumoniae_". Given almost 6 % of reads were "unclassified" and ~2 % of reads were assigned to very closely related taxa (in the _Klebsiella_ genus), this suggests the reads are from _Klebsiella pneumoniae_ with very little -if any- read contamination. 
    
    ```
     5.98	 108155	108155	U	0	unclassified
     94.02	1699669	0	C	1	
     94.02	1699669	1862	C1	131567	  cellular organisms
     93.91	1697788	2590	D	2	    Bacteria
     93.75	1694805	6312	P	1224	      Proteobacteria
     93.39	1688284	37464	C	1236	        Gammaproteobacteria
     91.31	1650648	35278	O	91347	          Enterobacterales
     89.31	1614639	43698	F	543	            Enterobacteriaceae
     86.40	1561902	22513	G	570	              Klebsiella
     **84.35	1524918	1524918	S	573	                Klebsiella pneumoniae**
      0.75	13596	13596	S	548	                Klebsiella aerogenes
      0.03	600	600	S	244366	                Klebsiella variicola
      0.01	253	253	S	571	                Klebsiella oxytoca
      0.00	17	17	S	1134687	                Klebsiella michiganensis
      0.00	3	0	G1	2608929	                unclassified Klebsiella
      0.00	3	3	S	1972757	                  Klebsiella sp. PO552
      0.00	2	2	S	1463165	                Klebsiella quasipneumoniae
      0.17	3035	129	G	590	              Salmonella
      0.15	2728	909	S	28901	                Salmonella enterica
      0.03	582	582	S1	9000010	                  Salmonella enterica subsp. IIa
      0.02	306	306	S1	59201	                  Salmonella enterica subsp. enterica
      0.01	230	230	S1	9000014	                  Salmonella enterica subsp. IIIa
      0.01	221	221	S1	9000015	                  Salmonella enterica subsp. IIIb
      0.01	136	136	S1	9000016	                  Salmonella enterica subsp. IX
      0.01	132	132	S1	9000011	                  Salmonella enterica subsp. IIb
      0.01	122	122	S1	59208	                  Salmonella enterica subsp. VII
      0.00	41	41	S1	59207	                  Salmonella enterica subsp. indica
      0.00	25	25	S1	9000017	                  Salmonella enterica subsp. X
      0.00	24	24	S1	9000009	                  Salmonella enterica subsp. VIII
      0.01	178	178	S	54736	                Salmonella bongori
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
