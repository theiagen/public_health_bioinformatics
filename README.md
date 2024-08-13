# Public Health Bioinformatics (PHB)

The Public Health Bioinformatics Bioinformatics repository contains workflows for genomic characterization, submission preparation, and genomic epidemiology of pathogens of public health concern.

## Introduction

**More information about the steps undertaken in these workflows is available via the [Theiagen Public Resources Documentation](https://theiagen.notion.site/Theiagen-Public-Health-Resources-a4bd134b0c5c4fe39870e21029a30566).**

Support for running these workflows can be sought by raising a [GitHub issue](https://github.com/theiagen/public_health_bioinformatics/issues/new) or by contacting Theiagen at support@theiagen.com.

These workflows are written in [WDL](https://github.com/openwdl/wdl), a language for specifying data processing workflows with a human-readable and writeable syntax. They have been developed by [Theiagen Genomics](https://theiagen.com/) to primarily run on the [Terra.bio](https://terra.bio/) platform but can be run locally or on an HPC system at the command-line with Cromwell or miniWDL.

## Purpose & Workflows

The PHB repository contains workflows for the characterization, genomic epidemiology, and sharing of pathogen genomes of public health concern. Workflows are available for viruses, bacteria, and fungi.

All workflows in the PHB repository end with `_PHB` in order to differentiate them from earlier versions and from the original tools they incorporate.

Briefly, the main *genomic characterization* workflows are split by pathogen type:

1. **Viral** (***TheiaCoV*** workflows)
2. **Bacterial** (***TheiaProk*** workflows)
3. **Fungal** (***TheiaEuk*** workflows)

Many more workflows are available, and are documented in detail in the [Theiagen Public Resources Documentation](https://theiagen.notion.site/Theiagen-Public-Health-Resources-a4bd134b0c5c4fe39870e21029a30566).

## On the Shoulder of Giants

The PHB repository would not be possible without its predecessors. We would like to acknowledge the following repositories, individuals, and contributors for their influence on the development of these workflows:

The PHB repository originated from collaborative work with Andrew Lang, PhD & his [Genomic Analysis WDL workflows](https://github.com/AndrewLangvt/genomic_analyses). The workflows and task development were influenced by The Broad's [Viral Pipes](https://github.com/broadinstitute/viral-pipelines) repository. The TheiaCoV workflows for viral genomic characterization were influenced by UPHL's [Cecret](https://github.com/UPHL-BioNGS/Cecret) & StaPH-B's [Monroe](https://staph-b.github.io/staphb_toolkit/workflow_docs/monroe/) workflows. The TheiaProk workflows for bacterial genomic characterization were influenced by Robert Petit's [bactopia](https://github.com/bactopia/bactopia). Most importantly, the PHB user community drove the development of these workflows and we are grateful for their feedback and contributions.

If you would like to provide feedback, please raise a [GitHub issue](https://github.com/theiagen/public_health_bioinformatics/issues/new).

### Contributing to the PHB workflows

We warmly welcome contributions to this repository! Our style guide may be found [here](https://theiagen.notion.site/Style-Guide-WDL-Workflow-Development-51b66a47dde54c798f35d673fff80249) for convenience of formatting.

You can expect a careful review of every PR and feedback as needed before merging, just like we do for PRs submitted by the Theiagen team. Our PR template can help prepare you for the review process. As always, reach out with any questions! We love recieving feedback and contributions from the community. When your PR is merged, we'll add your name to the contributors list below!

## Authorship & Responsibility

### Authorship

(Ordered by contribution [# of lines changed] as of 2024-03-28)

* **Sage Wright** ([@sage-wright](https://github.com/sage-wright)) - Conceptualization, Software, Validation, Supervision
* **Inês Mendes** ([@cimendes](https://github.com/cimendes)) - Software, Validation
* **Curtis Kapsak** ([@kapsakcj](https://github.com/kapsakcj)) - Conceptualization, Software, Validation
* **Frank Ambrosio** ([@frankambrosio3](https://github.com/frankambrosio3)) - Conceptualization, Software, Validation
* **James Otieno** ([@jrotieno](https://github.com/jrotieno)) - Software, Validation
* **Michelle Scribner** ([@michellescribner](https://github.com/michellescribner)) - Software, Validation
* **Kevin Libuit** ([@kevinlibuit](https://github.com/kevinlibuit)) - Conceptualization, Project Administration, Software, Validation, Supervision
* **Robert Petit** ([@rpetit3](https://github.com/rpetit3)) - Software, Validation
* **Emma Doughty** ([@emmadoughty](https://github.com/emmadoughty)) - Software, Validation
* **Andrew Page** ([@andrewjpage](https://github.com/andrewjpage)) - Project Administration, Software, Supervision
* **Kelsey Kropp** ([@kelseykropp](https://github.com/kelseykropp)) - Validation
* **Emily Smith** ([@emily-smith1](https://github.com/emily-smith1)) - Validation
* **Joel Sevinsky** ([@sevinsky](https://github.com/sevinsky)) - Conceptualization, Project Administration, Supervision

### External Contributors

We would like to gratefully acknowledge the following individuals from the public health community for their contributions to the PHB repository:

* **Ash O'Farrel** ([@aofarrel](https://github.com/aofarrel))
* **Holly Halstead** ([@HNHalstead](https://github.com/HNHalstead))
* **Sam Baird** ([@sam-baird](https://github.com/sam-baird))

### Maintaining PHB Pipelines

Theiagen Genomics has committed to maintaining these workflows for the forseeable future. These workflows are written using a standard workflow language (WDL) and uses Docker images based on the [StaPHB-B Docker Builds](https://github.com/StaPH-B/docker-builds). New versions that include bug fixes and additional features are released on a quarterly bases, with urgent bug fixes released as needed. Each version is accompanied by detailed release notes to lower the barrier of pipeline upkeep from the public health community at large.

### Point of Contact

If you have any questions or concerns, please raise a [GitHub issue](https://github.com/theiagen/public_health_bioinformatics/issues/new) or email Theiagen's general support at support@theiagen.com.

### Conflict of Interest

The authors declare no conflict of interest.

## Citation

Please cite this paper if publishing work using any workflows:

> Libuit, Kevin G., Emma L. Doughty, James R. Otieno, Frank Ambrosio, Curtis J. Kapsak, Emily A. Smith, Sage M. Wright, et al. 2023. “Accelerating Bioinformatics Implementation in Public Health.” Microbial Genomics 9 (7). https://doi.org/10.1099/mgen.0.001051.

Alternatively, please cite this paper if using the TheiaEuk workflow:

> Ambrosio, Frank, Michelle Scribner, Sage Wright, James Otieno, Emma Doughty, Andrew Gorzalski, Danielle Siao, et al. 2023. “TheiaEuk: A Species-Agnostic Bioinformatics Workflow for Fungal Genomic Characterization.” Frontiers in Public Health 11. https://doi.org/10.3389/fpubh.2023.1198213.
