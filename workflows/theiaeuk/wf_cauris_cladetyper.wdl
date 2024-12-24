version 1.0

import "../../tasks/species_typing/candida/task_cauris_cladetyper.wdl" as gambit_cladetyper_task
import "../../tasks/task_versioning.wdl" as versioning

workflow cauris_cladetyper {
  input {
    File assembly_fasta
    String samplename
  }
  call gambit_cladetyper_task.cauris_cladetyper as gambit_cladetyper {
    input:
      assembly_fasta = assembly_fasta,
      samplename = samplename
  }
  call versioning.version_capture {
    input:
  }
  output {
    String cauris_cladetyper_wf_version = version_capture.phb_version
    String cauris_cladetyper_wf_analysis_date = version_capture.date
    String gambit_cladetyper_clade = gambit_cladetyper.gambit_cladetype
    String gambit_cladetyper_annotated_reference = gambit_cladetyper.annotated_reference
    String gambit_cladetyper_version = gambit_cladetyper.version
  }
}