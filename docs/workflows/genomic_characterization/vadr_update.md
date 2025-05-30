# VADR_Update

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Genomic Characterization](../../workflows_overview/workflows_type.md/#genomic-characterization) | [Viral](../../workflows_overview/workflows_kingdom.md/#viral) | PHB v3.0.1 | Yes | Sample-level |

## Vadr_Update_PHB

The VADR_Update workflow updates prior VADR assessments for each sample in line with the assessment criteria in an alternative docker image. This may be useful when samples have previously been subject to VADR alerts as updates to VADR assessment criteria may mean that the sample no longer raises concern about quality. The latest docker image SARS-CoV-2 for VADR can be found [here](https://theiagen.notion.site/Docker-Image-and-Reference-Materials-for-SARS-CoV-2-Genomic-Characterization-98328c61f5cb4f77975f512b55d09108).

Various models are available for many organisms. The following table provides an overview of the recommended container to be used and what options should be passed on to VADR.

| **Organism** | **docker** | **vadr_opts** | max_length |
| --- | --- | --- | --- |
| sars-cov-2 | "us-docker.pkg.dev/general-theiagen/staphb/vadr:1.6.3" | "--noseqnamemax --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn --out_allfasta" | 30000 |
| MPXV | "us-docker.pkg.dev/general-theiagen/staphb/vadr:1.6.3" | "--glsearch -s -r --nomisc --mkey mpxv --r_lowsimok --r_lowsimxd 100 --r_lowsimxl 2000 --alt_pass discontn,dupregin --out_allfasta --minimap2 --s_overhang 150" | 210000 |
| WNV | "us-docker.pkg.dev/general-theiagen/staphb/vadr:1.6.3" | "--mkey flavi --mdir /opt/vadr/vadr-models-flavi/ --nomisc --noprotid --out_allfasta" | 11000 |
| flu | "us-docker.pkg.dev/general-theiagen/staphb/vadr:1.6.3" | "--atgonly --xnocomp --nomisc --alt_fail extrant5,extrant3 --mkey flu" | 13500 |
| rsv_a | "us-docker.pkg.dev/general-theiagen/staphb/vadr:1.6.3" | "-r --mkey rsv --xnocomp" | 15500 |
| rsv_b | "us-docker.pkg.dev/general-theiagen/staphb/vadr:1.6.3" | "-r --mkey rsv --xnocomp" | 15500 |
| HAV | "us-docker.pkg.dev/general-theiagen/staphb/vadr:1.6.3-hav" | "-r -xnocomp -mkey hav.vadr" | 10500 |

### Inputs

Please note the default values are for SARS-CoV-2.

This workflow runs on the sample level.

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "VADR_Update"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "VADR_Update"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///
