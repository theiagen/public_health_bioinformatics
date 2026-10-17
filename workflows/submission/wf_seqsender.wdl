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
    Boolean prep_only = false
    File? resume_submission_tarball # seqsender_state_tarball from an earlier run; skips build_and_prep and submit, and only polls
    Float max_wait_hours = 48
    Float poll_interval_minutes = 10

    File? input_table # for command line testing only; skips the Terra download
  }
  call versioning.version_capture {
    input:
  }
  if (!defined(resume_submission_tarball)) {
    call seqsender_tasks.seqsender_build_and_prep {
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
        prep_only = prep_only,
        input_table = input_table
    }
    # the tarball is only written if at least one sample passed prep
    if (!prep_only && defined(seqsender_build_and_prep.tarball)) {
      call seqsender_tasks.seqsender_submit {
        input:
          submission_tarball = select_first([seqsender_build_and_prep.tarball]),
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
  File? state_tarball = if defined(resume_submission_tarball) then resume_submission_tarball else seqsender_submit.tarball
  if (defined(state_tarball)) {
    call seqsender_tasks.seqsender_poll_and_wait {
      input:
        state_tarball = select_first([state_tarball]),
        seqsender_key = seqsender_key,
        max_wait_hours = max_wait_hours,
        poll_interval_minutes = poll_interval_minutes
    }
  }
  output {
    # Version Capture
    String seqsender_wf_version = version_capture.phb_version
    String seqsender_wf_analysis_date = version_capture.date
    # Build and Prep - skipped when resuming
    File? seqsender_metadata = seqsender_build_and_prep.combined_metadata
    File? seqsender_fasta = seqsender_build_and_prep.combined_fasta
    File? seqsender_build_and_prep_tarball = seqsender_build_and_prep.tarball
    File? seqsender_validation_summary = seqsender_build_and_prep.validation_summary
    Int? seqsender_num_samples_passed = seqsender_build_and_prep.num_samples_passed
    Int? seqsender_num_samples_failed = seqsender_build_and_prep.num_samples_failed
    File? seqsender_build_and_prep_log = seqsender_build_and_prep.log
    String? seqsender_version = seqsender_build_and_prep.seqsender_version
    String? seqsender_docker = seqsender_build_and_prep.seqsender_docker
    # Submit - skipped when resuming or with prep_only
    String? seqsender_submit_status = seqsender_submit.status
    # Poll - use `seqsender_state_tarball` with `resume_submission_tarball`
    File? seqsender_poll_tarball = seqsender_poll_and_wait.tarball
    String? seqsender_submission_name = seqsender_poll_and_wait.submission_name
    String? seqsender_final_status = seqsender_poll_and_wait.final_status
    String? seqsender_database_statuses = seqsender_poll_and_wait.database_statuses
    File? seqsender_poll_log = seqsender_poll_and_wait.log
    File? seqsender_submission_log = seqsender_poll_and_wait.submission_log
    File? seqsender_submission_status_report = seqsender_poll_and_wait.submission_status_report
    File? seqsender_biosample_report = seqsender_poll_and_wait.biosample_xml_report
    File? seqsender_sra_report = seqsender_poll_and_wait.sra_xml_report
    File? seqsender_genbank_report = seqsender_poll_and_wait.genbank_xml_report
  }
}
