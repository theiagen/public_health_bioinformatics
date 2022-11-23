version 1.0

import "../tasks/taxon_id/task_gambit.wdl" as gambit
import "../tasks/task_versioning.wdl" as versioning

workflow gambit_query {
  input {
    File assembly_fasta
    String samplename
  }
  call gambit.gambit {
    input:
      assembly = assembly_fasta,
      samplename = samplename,
  }
  call versioning.version_capture {
    input:
  }
  output {
    String gambit_query_wf_version = version_capture.phbg_version
    String gambit_query_wf_analysis_date = version_capture.date
    #Taxon ID
    File gambit_report = gambit.gambit_report_file
    File gambit_closest_genomes = gambit.gambit_closest_genomes_file
    String gambit_predicted_taxon = gambit.gambit_predicted_taxon
    String gambit_predicted_taxon_rank = gambit.gambit_predicted_taxon_rank
    String gambit_version = gambit.gambit_version
    String gambit_db_version = gambit.gambit_db_version
    String gambit_docker = gambit.gambit_docker
  }
}