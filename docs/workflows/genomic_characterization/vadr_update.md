# VADR_Update

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Name": "[**VADR_Update**](../workflows/genomic_characterization/vadr_update.md)"}, columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level", "Dockstore"]) }}

## VADR_Update_PHB

VADR_Update_PHB is a standalone workflow dedicated to running VADR. By default, the workflow uses a slimmed-down docker image running VADR (v1.6.4), which requires models to be provided separately. The table below outlines the recommended models and VADR parameters for use in the workflow.

| **Organism** | **vadr_model_file** | **vadr_opts** | **max_length** |
| --- | --- | --- | --- |
| sars-cov-2 | `"gs://theiagen-public-resources-rp/reference_data/databases/vadr_models/vadr-models-sarscov2-1.3-2.tar.gz"` | `"--mkey sarscov2 --glsearch -s -r --nomisc --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn --noseqnamemax --out_allfasta"` | `30000` |
| MPXV | `"gs://theiagen-public-resources-rp/reference_data/databases/vadr_models/vadr-models-mpxv-1.4.2-1.tar.gz"` | `"--mkey mpxv --glsearch --minimap2 -s -r --nomisc --r_lowsimok --r_lowsimxd 100 --r_lowsimxl 2000 --alt_pass discontn,dupregin --s_overhang 150 --out_allfasta"` | `210000` |
| WNV | `"gs://theiagen-public-resources-rp/reference_data/databases/vadr_models/vadr-models-flavi-1.2-1.tar.gz"` | `"--mkey flavi --nomisc --noprotid --out_allfasta"` | `11000` |
| flu | `"gs://theiagen-public-resources-rp/reference_data/databases/vadr_models/vadr-models-flu-1.6.3-2.tar.gz"` | `"--mkey flu --atgonly --xnocomp --nomisc --alt_fail extrant5,extrant3"` | `13500` |
| rsv_a | `"gs://theiagen-public-resources-rp/reference_data/databases/vadr_models/vadr-models-rsv-1.5-2.tar.gz"` | `"--mkey rsv --xnocomp -r"` | `15500` |
| rsv_b | `"gs://theiagen-public-resources-rp/reference_data/databases/vadr_models/vadr-models-rsv-1.5-2.tar.gz"` | `"--mkey rsv --xnocomp -r"` | `15500` |
| measles | `"gs://theiagen-public-resources-rp/reference_data/databases/vadr_models/vadr-models-mev-1.02.tar.gz"` | `"--mkey mev -r --indefclass 0.01"` | `18000` |
| mumps | `"gs://theiagen-public-resources-rp/reference_data/databases/vadr_models/vadr-models-muv-1.01.tar.gz"` | `"--mkey muv -r --indefclass 0.025"` | `18000` |
| rubella | `"gs://theiagen-public-resources-rp/reference_data/databases/vadr_models/vadr-models-ruv-1.01.tar.gz"` | `"--mkey ruv -r"` | `10000` |

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
