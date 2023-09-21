version 1.0

import "../../tasks/utilities/task_theiacov_fasta_utilities.wdl" as assembly_fasta_utils
import "../../tasks/task_versioning.wdl" as versioning

workflow theiacov_fasta_prep {
  input {
    File assembly
    String seq_method
    String assembly_method
    String organism
    String flu_segment
  }
  call assembly_fasta_utils.prep_theiacov_fasta_metadata {
    input:
      assembly = assembly,
      seq_method = seq_method,
      assembly_method = assembly_method,
      organism = organism,
      flu_segment = flu_segment
  }
  call versioning.version_capture {
    input:
  }
  output {
    String theaicov_fasta_set_prep_version = version_capture.phb_version
    String theaicov_fasta_set_prep_analysis_date = version_capture.date
    File theaicov_fasta_set_metadata = prep_theiacov_fasta_metadata.theaicov_fasta_set_metadata
  }
}