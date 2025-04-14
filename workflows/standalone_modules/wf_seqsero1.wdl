version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/species_typing/salmonella/task_seqsero1.wdl" as seqsero1

workflow seqsero1 {
  input {
    String samplename
    File assembly_fasta
  }
  call versioning.version_capture {
    input:
  }
  call seqsero1.seqsero`_assembly {
    input:
      samplename = samplename,
      assembly_fasta = assembly_fasta
  }
  output {
    String seqsero1_analysis_date = version_capture.date
    String seqsero1_wf_version = version_capture.phb_version
    String seqsero1_report = seqsero1.seqsero_report
    String seqsero1_version = seqsero1.seqsero_version
    String seqsero1_predicted_antigenic_profile = seqsero1.seqsero_predicted_antigenic_profile
    String seqsero1_predicted_serotype = seqsero1.seqsero_predicted_serotype
  }

}