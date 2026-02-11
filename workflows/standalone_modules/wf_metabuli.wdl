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
    Boolean call_trim = false
    Boolean illumina = false
  }
  call versioning_task.version_capture {
    input:
  }
  # Determine seq_mode
  if (! defined(read2) && illumina) {
    Int se_mode = 1
  }
  if (defined(read2)) {
    Int pe_mode = 2
  }
  if (! defined(read2) && ! illumina) {
    Int ont_mode = 3
  }
  # implicitly infer Illumina if read2 is provided
  if (defined(read2) && ! illumina) {
    Boolean implicit_illumina = true
    Int implicit_pe_mode = 2
  }
  # Trim reads if requested
  if (call_trim && (illumina || implicit_illumina)) {
    # Trim Illumina
    call fastp_task.fastp {
        input:
          read1 = read1,
          read2 = read2,
          samplename = samplename,
          fastp_trim_adapters = true
    }
  }
  # Trim ONT
  if (call_trim && !(illumina || implicit_illumina)) {
    call porechop_task.porechop {
      input:
        read1 = read1,
        samplename = samplename
    }
  }
  if (defined(read2)) {
    File read2_input = select_first([fastp.read2_trimmed, read2])
  }
  call metabuli_task.metabuli {
    input:
      samplename = samplename,
      read1 = select_first([fastp.read1_trimmed, porechop.trimmed_reads, read1]),
      read2 = read2_input,
      seq_mode = select_first([implicit_pe_mode, se_mode, pe_mode, ont_mode])
  }
  output {
    # PHB Version Captures
    String metabuli_wf_version = version_capture.phb_version
    String metabuli_wf_analysis_date = version_capture.date
    # Read trimming
    File? fastp_read1_trimmed = fastp.read1_trimmed
    File? fastp_read2_trimmed = fastp.read2_trimmed
    String? fastp_version = fastp.fastp_version
    String? fastp_docker = fastp.fastp_docker
    File? fastp_stats_html = fastp.fastp_stats_html
    File? fastp_stats_json = fastp.fastp_stats_json
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
