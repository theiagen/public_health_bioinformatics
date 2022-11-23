version 1.0

import "../tasks/species_typing/task_cauris_cladetyper.wdl" as gambit_cladetyper
import "../tasks/task_versioning.wdl" as versioning

workflow theiacauris_pe {
  input {
    File assembly_fasta
    String samplename
  }
  call gambit_cladetyper.cauris_cladetyper as gambit_cladetyper_task {
    input:
      assembly_fasta = assembly_fasta,
      samplename = samplename
  }
  call versioning.version_capture{
    input:
  }
  output {
    String theiacauris_pe_wf_version = version_capture.phbg_version
    String theiacauris_pe_wf_analysis_date = version_capture.date
    String theiacauris_pe_clade_assignment = gambit_cladetyper_task.gambit_cladetype
    String theiacauris_pe_docker = gambit_cladetyper_task.gambit_cladetyper_docker_image
  }
}