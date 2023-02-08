version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/task_pub_repo_prep.wdl" as submission_prep

workflow mercury_batch {
  input {
    Array[File] genbank_assembly
    Array[File] genbank_modifier
    Array[File] gisaid_assembly
    Array[File] gisaid_metadata
    Array[File] sra_metadata
    Array[String] sra_reads
    Array[File] biosample_attributes
    Array[String] samplename
    Array[String] submission_id
    Array[String] vadr_num_alerts
    Int vadr_threshold = 0
    Int cpu = 4
    Int disk_size = 100
    Int memory = 8
    String? gcp_bucket
  }
  call versioning.version_capture{
    input:
  }
  call submission_prep.compile_assembly_n_meta as genbank_compile {
    input:
      single_submission_fasta = genbank_assembly,
      single_submission_meta = genbank_modifier,
      samplename = samplename,
      vadr_num_alerts = vadr_num_alerts,
      repository = "GenBank",
      file_ext = "tsv",
      vadr_threshold = vadr_threshold,
      submission_id = submission_id,
      date = version_capture.date,
      cpu = cpu,
      disk_size = disk_size,
      memory = memory
  }
  call submission_prep.compile_assembly_n_meta as gisaid_compile {
    input:
      single_submission_fasta = gisaid_assembly,
      single_submission_meta = gisaid_metadata,
      samplename = samplename,
      vadr_num_alerts = vadr_num_alerts,
      repository = "GISAID",
      file_ext = "csv",
      vadr_threshold = vadr_threshold,
      submission_id = submission_id,
      date = version_capture.date,
      cpu = cpu,
      disk_size = disk_size,
      memory = memory
  }
  call submission_prep.compile_biosamp_n_sra {
    input:
      single_submission_biosample_attirbutes = biosample_attributes,
      single_submission_sra_metadata = sra_metadata,
      single_submission_sra_reads = sra_reads,
      gcp_bucket = gcp_bucket,
      date = version_capture.date,
      cpu = cpu,
      disk_size = disk_size,
      memory = memory
  }
  output {
    # Version Capture
    String mercury_batch_version = version_capture.phb_version
    String mercury_batch_analysis_date = version_capture.date
    # GenBank Submission Files
    File? GenBank_modifier  = genbank_compile.upload_meta
    File? GenBank_assembly = genbank_compile.upload_fasta
    File GenBank_batched_samples = genbank_compile.batched_samples
    File GenBank_excluded_samples = genbank_compile.excluded_samples
    # GISAID Submission Files
    File? GISAID_metadata  = gisaid_compile.upload_meta
    File? GISAID_assembly = gisaid_compile.upload_fasta
    File GISAID_batched_samples = gisaid_compile.batched_samples
    File GISAID_excluded_samples = gisaid_compile.excluded_samples
    # BioSample and SRA Submission Files
    File BioSample_attributes = compile_biosamp_n_sra.biosample_attributes
    File SRA_metadata = compile_biosamp_n_sra.sra_metadata
    File? SRA_zipped_reads = compile_biosamp_n_sra.sra_zipped
    String? SRA_gcp_bucket = gcp_bucket
  }
}