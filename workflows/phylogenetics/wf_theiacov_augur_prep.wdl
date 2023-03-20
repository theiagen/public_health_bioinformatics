version 1.0

import "../../tasks/phylogenetic_inference/task_nextstrain_augur.wdl" as nextstrain
import "../../tasks/task_versioning.wdl" as versioning

workflow theiacov_augur_prep {
  input {
    File assembly
    String collection_date
    String iso_country
    String iso_state
    String iso_continent
    String? iso_county
    String pango_lineage
  }
  call nextstrain.prep_augur_metadata {
    input:
      assembly = assembly,
      collection_date = collection_date,
      iso_country = iso_country,
      iso_state = iso_state,
      iso_continent = iso_continent,
      iso_county = iso_county,
      pango_lineage = pango_lineage
  }
  call versioning.version_capture {
    input:
  }
  output {
    String theiacov_augur_run_version = version_capture.phb_version
    String theiacov_augur_run_analysis_date = version_capture.date
    File augur_metadata = prep_augur_metadata.augur_metadata
  }
}