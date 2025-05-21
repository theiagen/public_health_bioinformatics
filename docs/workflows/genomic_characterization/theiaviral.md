# TheiaViral Workflow Series

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filter_column="Name", filter_values="[**TheiaViral Workflow Series**](../workflows/genomic_characterization/theiaviral.md)", columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level"]) }}

## TheiaViral Workflows

**TheiaViral** workflows assemble, quality assess, and characterize viral genomes from diverse data sources, including metagenomic samples. TheiaViral workflows can generate consensus assemblies of recalcitrant viruses, including diverse or recombinant lineages, such as rabies virus and norovirus, through a three-step approach: 1) generating an intermediate *de novo* assembly from taxonomy-filtered reads, 2) selecting the best reference from a database of 270,000+ complete viral genomes using average nucleotide identity, and 3) producing a final consensus assembly through reference-based read mapping and variant calling. While some sequencing methods, such as tiled amplicon sequencing, are incompatible with *de novo* assembly, TheiaViral may still be able to select a high-quality reference genome and generate a consensus assembly using these data sources. Targeted viral characterization is currently ongoing.

??? question "What are the main differences between the TheiaViral and TheiaCov workflows?"

    <div class="grid cards" markdown>

    -   :material-database: **TheiaCov Workflows**

        ---

        * For amplicon-derived viral sequencing methods
        * Supports a limited number of [pathogens](../../workflows/genomic_characterization/theiacov.md/#supported-organisms)
        * Uses manually curated, static reference genomes


    -   :material-database: **TheiaViral Workflows**

        ---

        * Designed for a variety of sequencing methods
        * Supports relatively diverse and recombinant pathogens
        * Dynamically identifies the most similar reference genome for consensus assembly via an intermediate *de novo* assembly

    </div>

### Workflow Diagram

=== "Illumina_PE"

    !!! caption "TheiaViral_Illumina_PE Workflow Diagram"

        ![TheiaViral_Illumina_PE Workflow Diagram](../../assets/figures/TheiaViral_Illumina_PE.png)

=== "ONT"

    !!! caption "TheiaViral_ONT Workflow Diagram"

        ![TheiaViral_ONT Workflow Diagram](../../assets/figures/TheiaViral_ONT.png)

## TheiaViral Workflows for Different Input Types

<div class="grid cards " markdown>

-   <center> **TheiaViral_Illumina_PE** </center>

    ---

    !!! dna "Illumina_PE Input Read Data"

        The TheiaViral_Illumina_PE workflow inputs Illumina paired-end read data. Read file extensions should be `.fastq` or `.fq`, and can optionally include the `.gz` compression extension. Theiagen recommends compressing files with [gzip](https://www.gnu.org/software/gzip/) before Terra uploads to minimize data upload time and storage costs.

        Modifications to the optional parameter for `trim_minlen` may be required to appropriately trim reads shorter than 2 x 150 bp (i.e. generated using a 300-cycle sequencing kit), such as the 2 x 75bp reads generated using a 150-cycle sequencing kit.

-   <center> **TheiaViral_ONT** </center>

    ---

    !!! dna_blue "ONT Input Read Data"

        The TheiaViral_ONT workflow inputs base-called Oxford Nanopore Technology (ONT) read data. Read file extensions should be `.fastq` or `.fq`, and can optionally include the `.gz` compression extension. Theiagen recommends compressing files with [gzip](https://www.gnu.org/software/gzip/) before Terra uploads to minimize data upload time and storage costs.

        It is recommended to trim adapter sequencings via `dorado` basecalling prior to running TheiaViral_ONT, though `porechop` can optionally be called to trim adapters within the workflow.

        **The ONT sequencing kit and base-calling approach can produce substantial variability in the amount and quality of read data. Genome assemblies produced by the TheiaViral_ONT workflow must be quality assessed before reporting results.**

</div>

### Inputs

=== "TheiaViral_Illumina_PE"
    /// html | div[class="searchable-table"]

    {{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="TheiaViral_Illumina_PE", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"], indent=8) }}

    ///

=== "TheiaViral_ONT"
    /// html | div[class="searchable-table"]

    {{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="TheiaViral_ONT", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"], indent=8) }}

    ///


#### Input notes

<div class="grid cards " markdown>
- ??? warning "`read_extraction_rank`"
    By default, the `read_extraction_rank` parameter is set to "family", which indicates that reads will be extracted if they are classified as the taxonomic family of the input `taxon`, including all descendant taxa of the family. Read classification may not resolve to the rank of the input `taxon`, so these reads may be classified at higher ranks. For example, some *Lyssavirus rabies* (species) reads may only be resolved to *Lyssavirus* (genus), so they would not be extracted if the `read_extraction_rank` is set to "species". Setting the `read_extraction_rank` above the inputted `taxon`'s rank can therefore dramatically increase the number of reads recovered, at the potential cost of including other viruses. This likely is not a problem for scarcely represented lineages, e.g. a sample that is expected to include *Lyssavirus rabies* is unlikely to contain other viruses of the corresponding family, Rhabdoviridae, within the same sample. However, setting a `read_extraction_rank` far beyond the input `taxon` rank can be problematic when multiple representatives of the same viral family are included in similar abundance within the same sample. 
    To further refine the desired `read_extraction_rank`, please review the corresponding classification reports of the respective classification software (kraken2 for Illumina and Metabuli for ONT)
</div>

<div class="grid cards " markdown>
- ??? warning "`min_allele_freq`, `min_depth`, and `min_map_quality`"
    These parameters have a direct effect on the variants that will ultimately be reported in the consensus assembly. `min_allele_freq` determines the minimum proportion of an allelic variant to be reported in the consensus assembly. `min_depth` and `min_map_quality` affect how "N" is reported in the consensus, i.e. depth below `min_depth` is reported as "N" and reads with mapping quality below `min_map_quality` are not included in depth calculations.      
</div>

### All Tasks

=== "TheiaViral_Illumina_PE"

    ??? toggle "Versioning"

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/versioning_task.md", condition="theiaviral", indent=12) }}

        </div>

    ??? toggle "Taxonomic Identification"

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/ncbi_identify_task.md", condition="theiaviral", indent=12) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/ncbi_taxon_summary_task.md", condition="theiaviral", indent=12) }}

        </div>

    ??? toggle "Read Quality Control, Trimming, Filtering, Identification and Extraction"

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/read_qc_trim_illumina.md", condition="theiaviral", indent=12, replacements={": Read Quality Trimming, Adapter Removal, Quantification, and Identification" : ""}) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/rasusa_task.md", condition="theiaviral", indent=12) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/read_screen_task.md", condition="theiaviral", indent=12, replacements={'??? task "`screen`: Total Raw Read Quantification and Genome Size Estimation"' : '??? task "`clean_check_reads`"'}) }}

        </div>

    ??? toggle "*De novo* Assembly and Reference Selection"

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/spades_task.md", condition="theiaviral", indent=12) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/megahit_task.md", condition="theiaviral", indent=12) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/skani_task.md", condition="theiaviral", indent=12) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/ncbi_datasets_task.md", condition="theiaviral", indent=12) }}

        </div>

    ??? toggle "Reference Mapping"

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/bwa_task.md", condition="theiaviral", indent=12) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/assembly_metrics_task.md", condition="theiaviral", indent=12, replacements={'`assembly_metrics`' : '`read_mapping_stats`'}) }}

        </div>

    ??? toggle "Variant Calling and Consensus Generation"

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/ivar_variants_task.md", condition="theiaviral", indent=12) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/ivar_consensus_task.md", condition="theiaviral", indent=12) }}

        </div>

    ??? toggle "Assembly Evaluation and Consensus Quality Control"

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/quast_task.md", condition="theiaviral", indent=12, replacements={'??? task "`quast`: Assembly Quality Assessment"' : '??? task "`quast_denovo`"'}) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/checkv_task.md", condition="theiaviral", indent=12, replacements={'??? task "`checkv`"' : '??? task "`checkv_denovo` & `checkv_consensus`"'}) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/consensus_qc_task.md", condition="theiaviral", indent=12) }}

        </div>

=== "TheiaViral_ONT"

    ??? toggle "Versioning"

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/versioning_task.md", condition="theiaviral", indent=12) }}

        </div>

    ??? toggle "Taxonomic Identification"

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/ncbi_identify_task.md", condition="theiaviral", indent=12) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/ncbi_taxon_summary_task.md", condition="theiaviral", indent=12) }}

        </div>

    ??? toggle "Read Quality Control, Trimming, and Filtering"

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/nanoplot_task.md", condition="theiaviral", indent=12, replacements={'??? task "`nanoplot`"' : '??? task "`nanoplot_raw` & `nanoplot_clean`"'}) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/porechop_task.md", condition="theiaviral", indent=12) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/nanoq_task.md", condition="theiaviral", indent=12) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/ncbi_scrub_task.md", condition="theiaviral", indent=12, replacements={'??? task "`HRRT`: Human Host Sequence Removal"' : '??? task "`ncbi_scrub_se`"'}) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/rasusa_task.md", condition="theiaviral", indent=12) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/read_screen_task.md", condition="theiaviral", indent=12, replacements={'??? task "`screen`: Total Raw Read Quantification and Genome Size Estimation"' : '??? task "`clean_check_reads`"'}) }}

        </div>

    ??? toggle "Read Classification and Extraction"

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/metabuli_task.md", condition="theiaviral", indent=12) }}

        </div>

    ??? toggle "*De novo* Assembly and Reference Selection"

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/raven_task.md", condition="theiaviral", indent=12) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/flye_task.md", condition="theiaviral", indent=12, replacements={'`flye_read_type`' : '`read_type`'}) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/skani_task.md", condition="theiaviral", indent=12) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/ncbi_datasets_task.md", condition="theiaviral", indent=12) }}

        </div>

    ??? toggle "Reference Mapping"

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/minimap2_task.md", condition="long_read_flags", indent=12, replacements={'??? task "`minimap2`: Read Alignment Details"' : '??? task "`minimap2`"'}) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/parse_mapping_task.md", condition="theiaviral_sam_to_sorted_bam", indent=12) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/assembly_metrics_task.md", condition="theiaviral", indent=12, replacements={'`assembly_metrics`' : '`read_mapping_stats`'}) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/fasta_utilities_task.md", condition="theiaviral", indent=12) }}

        </div>

    ??? toggle "Variant Calling and Consensus Generation"

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/clair3_task.md", condition="theiaviral", indent=12) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/parse_mapping_task.md", condition="theiaviral_mask_low_coverage", indent=12) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/bcftools_consensus_task.md", condition="theiaviral", indent=12) }}

        </div>

    ??? toggle "Assembly Evaluation and Consensus Quality Control"

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/quast_task.md", condition="theiaviral", indent=12, replacements={'??? task "`quast`: Assembly Quality Assessment"' : '??? task "`quast_denovo`"'}) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/checkv_task.md", condition="theiaviral", indent=12, replacements={'??? task "`checkv`"' : '??? task "`checkv_denovo` & `checkv_consensus`"'}) }}

        </div>

        <div class="grid cards" markdown>

        -   {{ include_md("common_text/consensus_qc_task.md", condition="theiaviral", indent=12) }}

        </div>

### Outputs

=== "TheiaViral_Illumina_PE"
    /// html | div[class="searchable-table"]

    {{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filter_column="Workflow", filter_values="TheiaViral_Illumina_PE", columns=["Variable", "Type", "Description"], sort_by=["Variable"], indent=8) }}

    ///

=== "TheiaViral_ONT"
    /// html | div[class="searchable-table"]

    {{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filter_column="Workflow", filter_values="TheiaViral_ONT", columns=["Variable", "Type", "Description"], sort_by=["Variable"], indent=8) }}

    ///

### Acknowlegments

We would like to thank Danny Park at the BROAD institute and Jared Johnson at the Washington State Department of Public Health for correspondence during the development of TheiaViral. TheiaViral was built referencing [viral-assemble](https://github.com/broadinstitute/viral-assemble/), [VAPER](https://github.com/DOH-JDJ0303/vaper), and [Artic](https://github.com/artic-network/fieldbioinformatics).