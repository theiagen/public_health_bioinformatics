version 1.0

import "../../tasks/utilities/task_augur_utilities.wdl" as augur_utils
import "../../tasks/task_versioning.wdl" as versioning

workflow augur_prep {
  input {
    String samplename
    String collection_date
    String country
    String state
    String continent
    String organism = "sars-cov-2" # options: "flu" or "sars-cov-2"
    String? pango_lineage
    String? county
  }
  call augur_utils.prep_augur_metadata {
    input:
      samplename = samplename,
      collection_date = collection_date,
      country = country,
      state = state,
      continent = continent,
      organism = organism,
      pango_lineage = pango_lineage,
      county = county
  }
  call versioning.version_capture {
    input:
  }
  output {
    String augur_prep_phb_version = version_capture.phb_version
    String augur_prep_phb_analysis_date = version_capture.date
    File augur_metadata = prep_augur_metadata.augur_metadata
  }
}