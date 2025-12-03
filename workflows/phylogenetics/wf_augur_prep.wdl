version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/data_handling/task_augur_utilities.wdl" as augur_utils

workflow augur_prep {
  input {
    File assembly
    String? collection_date
    String? country
    String? division
    String? region
    String? pango_lineage
    String? clade
    String? location
  }
  call augur_utils.prep_augur_metadata {
    input:
      assembly = assembly,
      collection_date = collection_date,
      country = country,
      division = division,
      region = region,
      pango_lineage = pango_lineage,
      nextclade_clade = clade,
      location = location
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