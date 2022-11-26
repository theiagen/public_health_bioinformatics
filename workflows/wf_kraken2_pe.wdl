version 1.0

import "../tasks/taxon_id/task_kraken2.wdl" as kraken2
import "../tasks/task_versioning.wdl" as versioning

workflow kraken2_pe_wf {
  meta {
    description: "Classify paired-end reads using Kraken2"
  }

  input {
    String  samplename
    File    read1
    File    read2
    File    kraken2_db
  }
  call kraken2.kraken2_pe {
    input:
      samplename = samplename,
      read1 = read1,
      read2 = read2,
      kraken2_db = kraken2_db
  }
  call versioning.version_capture{
    input:
  }
  output {
    # PHBG Version Captures
    String kraken2_pe_wf_version = version_capture.phbg_version
    String kraken2_pe_wf_analysis_date = version_capture.date
    # Kraken2
    String kraken2_version = kraken2_pe.kraken2_version
    String kraken2_docker = kraken2_pe.kraken2_docker
    File kraken2_report = kraken2_pe.kraken2_report
    File kraken2_classified_report = kraken2_pe.kraken2_classified_report
    File kraken2_unclassified_read1 = kraken2_pe.kraken2_unclassified_read1
    File kraken2_unclassified_read2 = kraken2_pe.kraken2_unclassified_read2
    File kraken2_classified_read1 = kraken2_pe.kraken2_classified_read1
    File kraken2_classified_read2 = kraken2_pe.kraken2_classified_read2
  }
}
