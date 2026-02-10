version 1.0

import "../../tasks/task_versioning.wdl" as versioning_task
import "../../tasks/taxon_id/contamination/task_metabuli.wdl" as metabuli_task
import "../../tasks/quality_control/read_filtering/task_fastp.wdl" as fastp_task
import "../../tasks/quality_control/read_filtering/task_porechop.wdl" as porechop_task

workflow metabuli_wf {
  meta {
    description: "Classify ONT/Illumina paired-end reads using Metabuli"
  }
  input {
    String samplename
    File read1
    File? read2
    Boolean illumina = false
    Boolean trim = false
  }
  call versioning_task.version_capture {
    input:
  }
  # Trim reads if requested
  if (trim) {
    # Trim Illumina PE
    if (defined(read2)) {
      call fastp_task.fastp_pe {
        input:
          read1 = read1,
          read2 = read2,
          samplename = samplename
      }
    }
    # Trim Illumina SE
    if (illumina && !defined(read2)) {
      call fastp_task.fastp_se {
        input:
          read1 = read1,
          samplename = samplename
      }
    }
    # Trim ONT
    if (!illumina && !defined(read2)) {
      call porechop_task.porechop {
        input:
          read1 = read1,
          samplename = samplename
      }
    }
  }
  call metabuli.metabuli as metabuli {
    input:
      samplename = samplename,
      read1 = select_first([fastp_pe.read1_trimmed, fastp_se.read1_trimmed, porechop.trimmed_reads, read1]),
      read2 = select_first([fastp_pe.read2_trimmed, fastp_se.read2_trimmed, read2])
  }
  output {
    # PHB Version Captures
    String metabuli_wf_version = version_capture.phb_version
    String metabuli_wf_analysis_date = version_capture.date
    # Read trimming
    File? fastp_read1_trimmed = select_first([fastp_pe.read1_trimmed, fastp_se.read1_trimmed])
    File? fastp_read2_trimmed = fastp_pe.read2_trimmed
    String? fastp_version = select_first([fastp_pe.fastp_version, fastp_se.fastp_version])
    String? fastp_docker = select_first([fastp_pe.fastp_docker, fastp_se.fastp_docker])
    File? fastp_stats_html = select_first([fastp_pe.fastp_stats_html, fastp_se.fastp_stats_html])
    File? fastp_stats_json = select_first([fastp_pe.fastp_stats_json, fastp_se.fastp_stats_json])
    File? porechop_read1_trimmed = porechop.trimmed_reads
    String? porechop_version = porechop.porechop_version
    # Metabuli
    String metabuli_version = metabuli.metabuli_version
    String metabuli_docker = metabuli.metabuli_docker
    File metabuli_report = metabuli.metabuli_report
    File metabuli_classified_report = metabuli.metabuli_classified
    File metabuli_krona_report = metabuli.metabuli_krona_report
    File? metabuli_classified_read1 = metabuli.metabuli_read1_extract
    File? metabuli_classified_read2 = metabuli.metabuli_read2_extract
  }
}
