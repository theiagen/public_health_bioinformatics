# Lyve_SET

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filter_column="Name", filter_values="[**Lyve_SET**](../workflows/phylogenetic_construction/lyve_set.md)", columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level"]) }}

## Lyve_SET_PHB

The Lyve_SET WDL workflow runs the [Lyve-SET](https://github.com/lskatz/lyve-SET) pipeline developed by Lee Katz et al. for phylogenetic analysis of bacterial genomes using high quality single nucleotide polymorphisms (hqSNPs). The Lyve_SET workflow identifies SNPs amongst a set of samples by mapping sequencing reads to a reference genome, identifying high quality SNPs, and inferring phylogeny using RAxML.

### Lyve-SET Pipeline (from [Lyve-SET paper](https://www.frontiersin.org/articles/10.3389/fmicb.2017.00375/full))

!!! caption "Lyve-SET Workflow Diagram"
    ![Lyve-SET Workflow Diagram](../../assets/figures/Lyve_Set.png)

### Inputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "Lyve_SET"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Workflow Actions

The Lyve_SET WDL workflow is run using read data from a set of samples. The workflow will produce a pairwise SNP matrix for the sample set and a maximum likelihood phylogenetic tree. Details regarding the default implementation of Lyve_SET and optional modifications are listed below.

1. Read processing
    1. By default, the Lyve_SET WDL workflow will perform read cleaning using the CG-Pipeline "CGP". However, read cleaning can be turned off or performed using "BayesHammer" using the `read_cleaner` input variable.
2. Reference procurement
    1. By default, the Lyve_SET WDL workflow will **not** mask phages or cliffs in the reference genome. Cliffs refer to regions of the reference genome where read coverage rises or falls abruptly. Masking phages and cliffs is intended to remove low quality SNPs. Users can invoke phage and cliff masking by setting the `mask_cliffs` and `mask_phages` variables to "true".
3. SNP discovery
    1. The Lyve_SET WDL workflow uses the default read mapper and variant caller from the Lyve-SET pipeline  (`smalt` and `varscan`). Additional options for each are available using the `mapper` and `snpcaller` input variables.
    2. The workflow also uses the default parameters for variant calling from the Lyve-SET pipeline: the minimum percent consensus to call a base is 0.75 and minimum read depth is 10X. These parameters can be manually modified using the `min_alt_frac` and `min_coverage` input variables.
4. Phylogenetic analysis
    1. The Lyve_SET workflow will attempt to produce a multiple sequence alignment, SNP distance matrix, and phylogenetic tree. These actions can be skipped by indicating `nomsa` = true, `nomatrix` = true, or `notrees` = true, respectively.

### Outputs

For full descriptions of Lyve-SET pipeline outputs, we recommend consulting the Lyve-SET documentation: <https://github.com/lskatz/lyve-SET/blob/master/docs/OUTPUT.md>

The following output files are populated to the Terra data table. However, please note that certain files may not appear in the data table following a run for two main reasons:

1. The user instructed the workflow to skip an analysis step
    1. For example, if `notrees` = true, no tree file will appear
2. The workflow skipped an analysis step due to an issue with the input data
    1. For example, the workflow will not attempt to produce a phylogenetic tree if there are too few samples or if samples are too closely related

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "Lyve_SET"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///

In addition to these outputs, all of the files produced by the Lyve-SET pipeline are available in the task-level outputs, including intermediate files and individual bam and vcf files for each sample. These files can be accessed viewing the execution directory for the run.

## References

> **Lyve-SET** Katz LS, Griswold T, Williams-Newkirk AJ, Wagner D, Petkau A, et al. (2017) A Comparative Analysis of the Lyve-SET Phylogenomics Pipeline for Genomic Epidemiology of Foodborne Pathogens. Frontiers in Microbiology 8.
