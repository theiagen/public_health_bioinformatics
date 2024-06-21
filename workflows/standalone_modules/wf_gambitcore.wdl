version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/quality_control/advanced_metrics/task_gambittools.wdl" as gambittools_task

workflow gambitcore_wf {
  input {
    File assembly_fasta
    String samplename
  }
  call gambittools_task.gambitcore {
    input:
      assembly = assembly_fasta,
      samplename = samplename
  }
  call versioning.version_capture {
    input:
  }
  output {
    String gambitcore_wf_version = version_capture.phb_version
    String gambitcore_wf_analysis_date = version_capture.date
    # Gambitcore output files
    File gambitcore_report_file = gambitcore.gambitcore_report_file
    String gambitcore_species = gambitcore.gambitcore_species
    String gambitcore_completeness = gambitcore.gambitcore_completeness
    String gambitcore_kmers_ratio = gambitcore.gambitcore_kmers_ratio
    String gambitcore_closest_accession = gambitcore.gambitcore_closest_accession
    String gambitcore_closest_distance = gambitcore.gambitcore_closest_distance
    String gambitcore_assembly_kmers = gambitcore.gambitcore_assembly_kmers
    String gambitcore_species_kmers = gambitcore.gambitcore_species_kmers
    String gambitcore_species_std_kmers = gambitcore.gambitcore_species_std_kmers
    String gambitcore_assembly_qc = gambitcore.gambitcore_assembly_qc
    String gambitcore_db_version = gambitcore.gambitcore_db_version
    String gambitcore_docker = gambitcore.gambitcore_docker
  }
}