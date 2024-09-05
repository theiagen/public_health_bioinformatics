# Usher

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line compatibility** | **Workflow type** |
|---|---|---|---|---|
| [Phylogenetic Placement](../../workflows_overview/workflows-type.md/#phylogenetic-placement) | [Viral](../../workflows_overview/workflows-kingdom.md/#viral) | PHB v2.1.0 | Yes | Sample-level, set-level |

## Usher_PHB

[UShER](https://usher-wiki.readthedocs.io/en/latest/) (Ultrafast Sample Placement on Existing Trees) rapidly places new samples onto an existing phylogeny using maximum parsimony. This workflow uses the UCSC-maintained global trees for SARS-CoV-2, mpox, RSV-A, and RSV-B if those organisms are specified in the `organism` input field. However, UShER can be used on any organism as long as a mutation-annotated tree (MAT) is provided in protobuf format. Contact us if you need help generating your own mutation-annotated tree, or follow the instructions available on the UShER wiki [here](https://usher-wiki.readthedocs.io/en/latest/).

### Inputs

While this workflow is technically a set-level workflow, it works on the sample-level too. When run on the set-level, the samples are placed with respect to each other.

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default attribute** | **Status** |
|---|---|---|---|---|---|
| usher_workflow | **assembly_fasta** | Array[File] | The assembly files for the samples you want to place on the pre-existing; can either be a set of samples, an individual sample, or multiple individual samples |  | Required |
| usher_workflow | **organism** | String | What organism to run UShER on; the following organism have default global phylogenies and reference files provided: sars-cov-2, mpox, RSV-A, RSV-B.  |  | Required |
| usher_workflow | **tree_name** | String | The output prefix for the uncondensed tree output and the clades output. |  | Required |
| usher | **cpu** | Int | CPUs to be allocated to UShER | 8 | Optional |
| usher | **disk_size** | Int | Disk size, in GB, to be allocated to UShER | 200 | Optional |
| usher | **docker** | String | The docker imaged used to run UShER | us-docker.pkg.dev/general-theiagen/pathogengenomics/usher:0.6.2 | Optional |
| usher | **memory** | Int | Memory, in GB, to be allocated to UShER | 32 | Optional |
| usher | **mutation_annotated_tree_pb** | File | Required for organisms other than sars-cov-2, mpox, RSV-A or RSV-B. This is the mutation-annotated global phylogeny upon which your samples will be placed  |  | Optional, Required |
| usher | **reference_genome** | File | Required for organisms other than sars-cov-2, mpox, RSV-A or RSV-B. This is the reference genome used to determine your sequence's mutations to accurately place the sample on the phylogeny.  |  | Optional, Required |
| usher | **subtree_size** | Int | Indicates how many of the closest-related samples you want to show in a subtree; more subtrees are made if there is more sequence diversity in the set of input samples (multiple subtrees are only generated if this workflow is run on the set level). | 20 | Optional |
| version_capture | **docker** | String | The Docker image used to run the version_capture task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) |  | Optional |

### Outputs

| **Variable** | **Type** | **Description** |
|---|---|---|
| usher_clades | File | The clades predicted for the samples |
| usher_phb_analysis_date | String | The date the analysis was run |
| usher_phb_version | String | The version of PHB the workflow is from |
| usher_protobuf_version | String | The version of the mutation-annotated protobuf tree (what day and what samples are included, if a default organism was used; otherwise, says it was user-provided) |
| usher_subtree_mutations | Array[File] | An array of files showing the mutations at each internal node for the subtree |
| usher_subtrees | Array[File] | An array of subtrees where your samples have been placed |
| usher_uncondensed_tree | File | The entire global tree with your samples included (warning: may be a very large file if the organism is "sars-cov-2") |
| usher_version | String | The version of UShER used |
