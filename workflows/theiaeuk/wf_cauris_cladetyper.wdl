version 1.0

import "../../tasks/species_typing/candidozyma/task_cauris_cladetyper.wdl" as gambit_cladetyper
import "../../tasks/task_versioning.wdl" as versioning

workflow cauris_cladetyper {
  input {
    File assembly_fasta
    String samplename
  }
  call gambit_cladetyper.cauris_cladetyper as cladetyper {
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
    String cladetyper_clade = cladetyper.gambit_cladetype
    String cladetyper_annotated_reference = cladetyper.annotated_reference
    String cladetyper_gambit_version = cladetyper.gambit_version
    String cladetyper_docker_image = cladetyper.gambit_cladetyper_docker_image
  }
}