version 1.0

import "../../tasks/task_versioning.wdl" as versioning_task
import "../../tasks/taxon_id/contamination/task_metabuli.wdl" as metabuli_task
import "../../tasks/quality_control/read_filtering/task_fastp.wdl" as fastp_task
import "../../tasks/quality_control/read_filtering/task_fastplong.wdl" as fastplong_task
import "../../tasks/taxon_id/task_ete4_taxon_id.wdl" as identify_taxon_id_task

workflow metabuli_wf {
  meta {
    description: "Classify ONT/Illumina paired-end reads using Metabuli"
  }
  input {
    String samplename
    String? taxon
    File read1
    File? read2
    File metabuli_db
    Boolean call_trim = true
    Boolean? illumina
    Int metabuli_mem = 32
    Int metabuli_disk_size = 250
  }
  call versioning_task.version_capture {
    input:
  }
  if (defined(taxon)) {
    call identify_taxon_id_task.ete4_taxon_id as ete4_identify {
      input:
        taxon = select_first([taxon])
    }
  }
  # Determine seq_mode
  if (! defined(read2) && select_first([illumina, false])) {
    Int se_mode = 1
  }
  if (defined(read2)) {
    Int pe_mode = 2
    Boolean implicit_illumina = true
  }
  if (! defined(read2) && ! select_first([illumina, false])) {
    Int ont_mode = 3
  }
  # Trim reads if requested
  if (call_trim && (select_first([illumina, false]) || defined(implicit_illumina))) {
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
  if (call_trim && ! (select_first([illumina, false]) || defined(implicit_illumina))) {
    call fastplong_task.fastplong {
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
      taxon_id = ete4_identify.taxon_id,
      read1 = select_first([fastp.read1_trimmed, fastplong.read1_trimmed, read1]),
      read2 = read2_input,
      seq_mode = select_first([se_mode, pe_mode, ont_mode]),
      metabuli_db = metabuli_db,
      memory = metabuli_mem,
      disk_size = metabuli_disk_size
  }
  output {
    # PHB Version Captures
    String metabuli_wf_version = version_capture.phb_version
    String metabuli_wf_analysis_date = version_capture.date
    # Taxon ID
    String? ncbi_taxon_id = ete4_identify.taxon_id
    String? ncbi_taxon_name = ete4_identify.taxon_name
    String? ncbi_read_extraction_rank = ete4_identify.taxon_rank
    String? ete4_version = ete4_identify.ete4_version
    String? ete4_docker = ete4_identify.ete4_docker
    # Read trimming
    File? fastp_read1_trimmed = fastp.read1_trimmed
    File? fastp_read2_trimmed = fastp.read2_trimmed
    String? fastp_version = fastp.fastp_version
    String? fastp_docker = fastp.fastp_docker
    File? fastp_html_report = fastp.fastp_stats_html
    File? fastp_json_report = fastp.fastp_stats_json
    File? fastplong_read1_trimmed = fastplong.read1_trimmed
    File? fastplong_html_report = fastplong.fastplong_stats_html
    File? fastplong_json_report = fastplong.fastplong_stats_json
    String? fastplong_version = fastplong.fastplong_version
    String? fastplong_docker = fastplong.fastplong_docker
    # Metabuli
    String metabuli_version = metabuli.metabuli_version
    String metabuli_docker = metabuli.metabuli_docker
    File metabuli_report = metabuli.metabuli_report
    File metabuli_classified_report = metabuli.metabuli_classified
    File metabuli_krona_report = metabuli.metabuli_krona_report
    File? metabuli_classified_read1 = metabuli.metabuli_read1_extract
    File? metabuli_classified_read2 = metabuli.metabuli_read2_extract
    String metabuli_status = metabuli.metabuli_status
  }
}
