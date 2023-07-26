# Public Health Bioinformatics (PHB)
Bioinformatics workflows for characterization, epidemiology and sharing of pathogen genomes.

**More information about the steps undertaken in these workflows is available via the [Theiagen Public Resources Documentation](https://theiagen.notion.site/Theiagen-Public-Health-Resources-a4bd134b0c5c4fe39870e21029a30566).**

Support for running these workflows can be sought by raising a [GitHub issue](https://github.com/theiagen/public_health_bioinformatics/issues/new) or by contacting Theiagen at support@theiagen.com.

These workflows are written in [WDL](https://github.com/openwdl/wdl), a language for specifying data processing workflows with a human-readable and writeable syntax. They have been developed by [Theiagen Genomics](https://theiagen.com/) to primarily run on the [Terra.bio](https://terra.bio/) platform but can be run locally or on an HPC system at the command-line with Cromwell or miniWDL.

### Contributors & Influence
* Based on collaborative work with Andrew Lang, PhD & his [Genomic Analysis WDL workflows](https://github.com/AndrewLangvt/genomic_analyses)
* Workflows and task development influenced by The Broad's [Viral Pipes](https://github.com/broadinstitute/viral-pipelines)
* TheiaCoV workflows for viral genomic characterization influenced by UPHL's [Cecret](https://github.com/UPHL-BioNGS/Cecret) & StaPH-B's [Monroe](https://staph-b.github.io/staphb_toolkit/workflow_docs/monroe/)
* TheiaProk workflows for bacterial genomic characterization influenced by Robert Petit's [bactopia](https://github.com/bactopia/bactopia)
* The PHB workflow user community. To provide feedback, please raise a [GitHub issue](https://github.com/theiagen/public_health_vioinformatics/issues/new).

### Contributing to the PHB workflows
Contributions to the workflows contained in this repository are warmly welcomed. Our style guide may be found [here](https://theiagen.notion.site/Style-Guide-WDL-Workflow-Development-bb456f34322d4f4db699d4029050481c) for convenience of formatting.

### Citation

Please cite this paper if publishing work using any workflows:

Libuit, Kevin G., Emma L. Doughty, James R. Otieno, Frank Ambrosio, Curtis J. Kapsak, Emily A. Smith, Sage M. Wright, et al. 2023. “Accelerating Bioinformatics Implementation in Public Health.” Microbial Genomics 9 (7). https://doi.org/10.1099/mgen.0.001051.

Alternatively, please cite this paper if using the TheiaEuk workflow:

Ambrosio, Frank, Michelle Scribner, Sage Wright, James Otieno, Emma Doughty, Andrew Gorzalski, Danielle Siao, et al. 2023. “TheiaEuk: A Species-Agnostic Bioinformatics Workflow for Fungal Genomic Characterization.” Frontiers in Public Health 11. https://doi.org/10.3389/fpubh.2023.1198213.