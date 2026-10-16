version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/submission/task_seqsender.wdl" as seqsender_tasks

workflow seqsender {
  input {
    String project_name
    String workspace_name
    String table_name
    File seqsender_config # NCBI password already encrypted with `seqsender.py load_credentials`
    String seqsender_key # the key that encrypted the password in seqsender_config
    Array[String] sample_ids

    String bioproject
    Array[File]? assembly_fastas # required for GENBANK
    Array[File]? read1_fastqs # required for SRA
    Array[File]? read2_fastqs # optional for SRA (single-end)
    String? sra_transfer_gcp_bucket # NCBI-readable bucket the passing samples' reads are copied to

    String organism = "FLU" # options: "FLU", "COV", "OTHER"
    Boolean submit_to_biosample = true
    Boolean submit_to_sra = true
    Boolean submit_to_genbank = true

    Boolean submit_to_production = false
    Boolean prep_only = false

    File? resume_submission_tarball # seqsender_state_tarball from an earlier run; skips build_and_prep and submit, and only polls

    File? input_table # for command line testing only; skips the Terra download

  }
  call versioning.version_capture {
    input:
  }
  if (!defined(resume_submission_tarball)) {
    call seqsender_tasks.seqsender_build_and_prep {
      input:
        project_name = project_name,
        workspace_name = workspace_name,
        table_name = table_name,
        sample_ids = sample_ids,
        seqsender_config = seqsender_config,
        bioproject = bioproject,
        sra_transfer_gcp_bucket = sra_transfer_gcp_bucket,
        organism = organism,
        submit_to_biosample = submit_to_biosample,
        submit_to_sra = submit_to_sra,
        submit_to_genbank = submit_to_genbank,
        assembly_fastas = assembly_fastas,
        read1_fastqs = read1_fastqs,
        read2_fastqs = read2_fastqs,
        input_table = input_table
    }
    # the tarball is only written if at least one sample passed prep
    if (!prep_only && defined(seqsender_build_and_prep.tarball)) {
      call seqsender_tasks.seqsender_submit {
        input:
          submission_tarball = select_first([seqsender_build_and_prep.tarball]),
          seqsender_key = seqsender_key,
          submit_to_production = submit_to_production,
          organism = organism,
          submit_to_biosample = submit_to_biosample,
          submit_to_sra = submit_to_sra,
          submit_to_genbank = submit_to_genbank,
      }
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
    File? seqsender_submit_log = seqsender_submit.log
  }
}
