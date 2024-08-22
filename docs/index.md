---
title: Overview
---

## Documentation and Resources

Theiagen’s approach to genomic analysis in public health typically uses the [Terra](https://terra.bio/) platform to run workflows that undertake bioinformatic analysis, then uses other platforms for visualization of the resulting data. Theiagen’s **Public Health Bioinformatics (PHB)** repository contains workflows for the characterization, genomic epidemiology, and sharing of pathogen genomes of public health concern. Workflows are available for viruses, bacteria, and fungi.

!!! dna "Our Open Source Philosophy"
    PHB source code is publicly available on [GitHub](https://github.com/theiagen/public_health_bioinformatics)!

    All workflows can be imported directly to [Terra](https://terra.bio/) via the [**Dockstore PHB collection**](https://dockstore.org/organizations/Theiagen/collections/public-health-bioinformatics)!

**When undertaking genomic analysis using Terra and other data visualization platforms, it is essential to consider the necessary and appropriate workflows and resources for your analysis. To help you make these choices, take a look at the relationship between the most commonly used Theiagen workflows.**

!!! caption "Analysis Approaches for Genomic Data"
    ![The relationship between the various PHB workflows](assets/Workflow_Relationships.png#only-light)
    ![The relationship between the various PHB workflows](assets/Workflow_Relationships_dark.png#only-dark)

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
