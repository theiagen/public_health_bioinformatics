version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/submission/task_seqsender.wdl" as seqsender_tasks

workflow seqsender {
  input {
    String terra_project_name
    String terra_workspace_name
    String terra_table_name
    File seqsender_config # NCBI password already encrypted with `seqsender.py load_credentials`
    String seqsender_key # the key that encrypted the password in seqsender_config
    Array[String] sample_ids

    String bioproject
    Array[File]? assembly_fastas # required for GENBANK
    Array[File]? read1_fastqs # required for SRA
    Array[File]? read2_fastqs # optional for SRA (single-end)
    String? sra_transfer_gcp_bucket # NCBI-readable bucket the passing samples' reads are copied to

    String seqsender_organism_tag # options: "FLU", "COV", "OTHER"
    Boolean submit_to_biosample
    Boolean submit_to_sra
    Boolean submit_to_genbank

    Boolean submit_to_production = false
    Boolean validation_only = false
    File? resume_checkpoint_tarball # checkpoint_tarball from an earlier run; skips build_and_prep and submit, and only polls
    Float max_wait_hours = 48
    Float poll_interval_minutes = 10

    File? input_table # for command line testing only; skips the Terra download
  }
  call versioning.version_capture {
    input:
  }
  if (!defined(resume_checkpoint_tarball)) {
    call seqsender_tasks.seqsender_build_and_prep as validation {
      input:
        terra_project_name = terra_project_name,
        terra_workspace_name = terra_workspace_name,
        terra_table_name = terra_table_name,
        seqsender_config = seqsender_config,
        sample_ids = sample_ids,
        bioproject = bioproject,
        assembly_fastas = assembly_fastas,
        read1_fastqs = read1_fastqs,
        read2_fastqs = read2_fastqs,
        sra_transfer_gcp_bucket = sra_transfer_gcp_bucket,
        seqsender_organism_tag = seqsender_organism_tag,
        submit_to_biosample = submit_to_biosample,
        submit_to_sra = submit_to_sra,
        submit_to_genbank = submit_to_genbank,
        validation_only = validation_only,
        input_table = input_table
    }
    # the tarball is only written if at least one sample passed prep
    if (!validation_only && defined(validation.tarball)) {
      call seqsender_tasks.seqsender_submit as submission {
        input:
          submission_tarball = select_first([validation.tarball]),
          seqsender_key = seqsender_key,
          submit_to_production = submit_to_production,
          seqsender_organism_tag = seqsender_organism_tag,
          submit_to_biosample = submit_to_biosample,
          submit_to_sra = submit_to_sra,
          submit_to_genbank = submit_to_genbank,
      }
    }
  }
  # poll a new submission, or resume polling an earlier one
  File? latest_tarball = if defined(resume_checkpoint_tarball) then resume_checkpoint_tarball else submission.tarball
  if (defined(latest_tarball)) {
    call seqsender_tasks.seqsender_poll_and_wait as polling {
      input:
        checkpoint_tarball = select_first([latest_tarball]),
        seqsender_key = seqsender_key,
        max_wait_hours = max_wait_hours,
        poll_interval_minutes = poll_interval_minutes
    }
  }
  output {
    # Version Capture
    String seqsender_wf_version = version_capture.phb_version
    String seqsender_wf_analysis_date = version_capture.date
    String? seqsender_version = select_first([validation.seqsender_version, polling.seqsender_version])
    String? seqsender_docker = select_first([validation.seqsender_docker, polling.seqsender_docker])
    # Submitted files - from the checkpoint when polling ran (incl. resuming), otherwise from validation; both are null when no sample passed
    File? submission_metadata = if defined(polling.submission_metadata) then polling.submission_metadata else validation.combined_metadata
    File? submission_fasta = if defined(polling.submission_fasta) then polling.submission_fasta else validation.combined_fasta
    # Build and Prep - skipped when resuming
    Int? validation_samples_passed = validation.num_samples_passed
    Int? validation_samples_failed = validation.num_samples_failed
    File? validation_summary = validation.validation_summary
    # Poll - use `checkpoint_tarball` with `resume_checkpoint_tarball`
    File? checkpoint_tarball = polling.tarball
    String? seqsender_final_status = polling.final_status
    File? seqsender_database_status = polling.submission_log_csv
    File? seqsender_submission_status = polling.submission_status_report
    File? seqsender_report_xml_biosample = polling.report_xml_biosample
    File? seqsender_report_xml_sra = polling.report_xml_sra
    File? seqsender_report_xml_genbank = polling.report_xml_genbank
    # Logs from each task
    File? log_validation = validation.log
    File? log_submission = submission.log
    File? log_polling = polling.log
  }
}
