version 1.0

import "../../tasks/taxon_id/freyja_utilities/task_freyja_update.wdl" as update
import "../../tasks/taxon_id/freyja_utilities/task_freyja_transfer.wdl" as transfer

workflow freyja_update {
  input {
    String gcp_uri
  }
  call update.freyja_update_refs {
    input:
  }
  call transfer.transfer_files {
    input:
      updated_barcodes = freyja_update_refs.updated_barcodes,
      updated_lineages = freyja_update_refs.updated_lineages,
      update_log = freyja_update_refs.update_log,
      gcp_uri = gcp_uri
  }
  output {
  }
}
