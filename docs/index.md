---
title: Overview
---

## Documentation and Resources

Theiagen's approach to genomic analysis in public health typically uses the [Terra](https://terra.bio/) platform to run workflows that undertake bioinformatic analysis, then uses other platforms for visualization of the resulting data. Theiagen's **Public Health Bioinformatics (PHB)** repository contains workflows for the characterization, genomic epidemiology, and sharing of pathogen genomes of public health concern. Workflows are available for viruses, bacteria, and fungi.

!!! dna "Our Open Source Philosophy"
    PHB source code is publicly available on [GitHub](https://github.com/theiagen/public_health_bioinformatics)!

    All workflows can be imported directly to [Terra](https://terra.bio/) via the [**Dockstore PHB collection**](https://dockstore.org/organizations/Theiagen/collections/public-health-bioinformatics)!

**When undertaking genomic analysis using Terra and other data visualization platforms, it is essential to consider the necessary and appropriate workflows and resources for your analysis. To help you make these choices, take a look at the relationship between the most commonly used Theiagen workflows.**

!!! caption "Analysis Approaches for Genomic Data"
    ![The relationship between the various PHB workflows](assets/figures/Workflow_Relationships.png#only-light){data-description="This diagram shows the Theiagen workflows (green boxes) available for analysis of genomic data in public health and the workflows that may be used consecutively (arrows). The blue boxes describe the major functions that these workflows undertake. The yellow boxes show functions that may be undertaken independently of workflows on Terra."}
    ![The relationship between the various PHB workflows](assets/figures/Workflow_Relationships_dark.png#only-dark){data-description="This diagram shows the Theiagen workflows (green boxes) available for analysis of genomic data in public health and the workflows that may be used consecutively (arrows). The blue boxes describe the major functions that these workflows undertake. The yellow boxes show functions that may be undertaken independently of workflows on Terra."}

    This diagram shows the Theiagen workflows (green boxes) available for analysis of genomic data in public health and the workflows that may be used consecutively (arrows). The blue boxes describe the major functions that these workflows undertake. The yellow boxes show functions that may be undertaken independently of workflows on Terra.

### PHB development is a cycle

We continuously work to improve our codebase and usability of our workflows by the public health community, so changes from version to version are expected.  This documentation page reflects the state of the workflow at the version stated in the title.

!!! dna "What's new?"
    You can see the changes since PHB v2.0.1 [**here**](https://www.notion.so/Public-Health-Bioinformatics-v2-1-0-Minor-Release-Notes-e4eb83259b744ca591fc5a6c4d53f977?pvs=21)!

Our code is open-source and available under [GNU Affero General Public License v3.0](https://github.com/theiagen/public_health_viral_genomics/blob/main/LICENSE)

----

### Contributing to the PHB Repository

If you would like to submit suggested code changes to our workflows, you may add or modify the WDL files and submit pull requests to the [PHB GitHub](https://github.com/theiagen/public_health_bioinformatics) repository.

### About Theiagen

Theiagen develops bioinformatics solutions for public health labs, and then trains and supports scientists to use these. If you would like to work with Theiagen, please [get in contact](https://theiagen.com/team-up-with-theiagen/).

------

The PHB repository contains workflows for the characterization, genomic epidemiology, and sharing of pathogen genomes of public health concern. Workflows are available for viruses, bacteria, and fungi.

All workflows in the PHB repository end with _PHB in order to differentiate them from earlier versions and from the original tools they incorporate.

Briefly, the main genomic characterization workflows are split by pathogen type:

Viral (TheiaCoV workflows)
Bacterial (TheiaProk workflows)
Fungal (TheiaEuk workflows)
Many more workflows are available, and are documented in detail in the Theiagen Public Resources Documentation.

On the Shoulder of Giants
The PHB repository would not be possible without its predecessors. We would like to acknowledge the following repositories, individuals, and contributors for their influence on the development of these workflows:

The PHB repository originated from collaborative work with Andrew Lang, PhD & his Genomic Analysis WDL workflows. The workflows and task development were influenced by The Broad's Viral Pipes repository. The TheiaCoV workflows for viral genomic characterization were influenced by UPHL's Cecret & StaPH-B's Monroe workflows. The TheiaProk workflows for bacterial genomic characterization were influenced by Robert Petit's bactopia. Most importantly, the PHB user community drove the development of these workflows and we are grateful for their feedback and contributions.

If you would like to provide feedback, please raise a GitHub issue.

Contributing to the PHB workflows
We warmly welcome contributions to this repository! Our style guide may be found here for convenience of formatting.

You can expect a careful review of every PR and feedback as needed before merging, just like we do for PRs submitted by the Theiagen team. Our PR template can help prepare you for the review process. As always, reach out with any questions! We love recieving feedback and contributions from the community. When your PR is merged, we'll add your name to the contributors list below!

Authorship & Responsibility
Authorship
(Ordered by contribution [# of lines changed] as of 2024-08-01)

Sage Wright (@sage-wright) - Conceptualization, Software, Validation, Supervision
Inês Mendes (@cimendes) - Software, Validation
Curtis Kapsak (@kapsakcj) - Conceptualization, Software, Validation
James Otieno (@jrotieno) - Software, Validation
Frank Ambrosio (@frankambrosio3) - Conceptualization, Software, Validation
Michelle Scribner (@michellescribner) - Software, Validation
Kevin Libuit (@kevinlibuit) - Conceptualization, Project Administration, Software, Validation, Supervision
Emma Doughty (@emmadoughty) - Software, Validation
Andrew Page (@andrewjpage) - Project Administration, Software, Supervision
Andrew Lang (@AndrewLangVt) - Software, Supervision
Kelsey Kropp (@kelseykropp) - Validation
Emily Smith (@emily-smith1) - Validation
Joel Sevinsky (@sevinsky) - Conceptualization, Project Administration, Supervision
External Contributors
We would like to gratefully acknowledge the following individuals from the public health community for their contributions to the PHB repository:

Robert Petit (@rpetit3)
Ash O'Farrel (@aofarrel)
Sam Baird (@sam-baird)
Holly Halstead (@HNHalstead)
Maintaining PHB Pipelines
Theiagen Genomics has committed to maintaining these workflows for the forseeable future. These workflows are written using a standard workflow language (WDL) and uses Docker images based on the StaPHB-B Docker Builds. New versions that include bug fixes and additional features are released on a quarterly bases, with urgent bug fixes released as needed. Each version is accompanied by detailed release notes to lower the barrier of pipeline upkeep from the public health community at large.

Point of Contact
If you have any questions or concerns, please raise a GitHub issue or email Theiagen's general support at support@theiagen.com.

Conflict of Interest
The authors declare no conflict of interest.

Citation
Please cite this paper if publishing work using any workflows:

Libuit, Kevin G., Emma L. Doughty, James R. Otieno, Frank Ambrosio, Curtis J. Kapsak, Emily A. Smith, Sage M. Wright, et al. 2023. “Accelerating Bioinformatics Implementation in Public Health.” Microbial Genomics 9 (7). https://doi.org/10.1099/mgen.0.001051.

Alternatively, please cite this paper if using the TheiaEuk workflow:

Ambrosio, Frank, Michelle Scribner, Sage Wright, James Otieno, Emma Doughty, Andrew Gorzalski, Danielle Siao, et al. 2023. “TheiaEuk: A Species-Agnostic Bioinformatics Workflow for Fungal Genomic Characterization.” Frontiers in Public Health 11. https://doi.org/10.3389/fpubh.2023.1198213.