version 1.0

import "../../../tasks/phylogenetic_inference/task_referenceseeker.wdl" as referenceseeker_task
import "../../../tasks/utilities/task_ncbi_datasets.wdl" as ncbi_datasets_task
import "../../../tasks/task_versioning.wdl" as versioning

# input is user-defined bacterial assembly_fasta file
# reference seeker takes in one sample
# get referenceseeker genome accession
# use ncbi datasets task to download reference genome fasta file

workflow reference_fetch {
  input {
    String samplename
    File? assembly_fasta
    String? ncbi_accession
  }
  # if user does not provide reference genome fasta, determine one for the user by running referenceseeker and ncbi datasets to acquire one
  if(! defined(ncbi_accession) && defined(assembly_fasta)){
  call referenceseeker_task.referenceseeker {
    input:
      assembly_fasta = select_first([assembly_fasta,""]), # to avoid error with required task input
      samplename = samplename
  }
  }
  call ncbi_datasets_task.ncbi_datasets_download_genome_accession {
    input:
      ncbi_accession = select_first([ncbi_accession,referenceseeker.referenceseeker_top_hit_ncbi_accession])
  }
  call versioning.version_capture {
    input:
  }
  output {
    # version capture
    String reference_fetch_version = version_capture.phb_version
    String reference_fetch_analysis_date = version_capture.date
    # referenceseeker outputs
    String? reference_fetch_referenceseeker_top_hit_ncbi_accession = referenceseeker.referenceseeker_top_hit_ncbi_accession
    String? reference_fetch_referenceseeker_version = referenceseeker.referenceseeker_version
    File? reference_fetch_referenceseeker_tsv = referenceseeker.referenceseeker_tsv
    String? reference_fetch_referenceseeker_docker = referenceseeker.referenceseeker_docker
    String? reference_fetch_referenceseeker_database = referenceseeker.referenceseeker_database
    # ncbi datasets outputs
    File reference_fetch_ncbi_datasets_assembly_data_report_json = ncbi_datasets_download_genome_accession.ncbi_datasets_assembly_data_report_json
    File reference_fetch_ncbi_datasets_assembly_fasta = ncbi_datasets_download_genome_accession.ncbi_datasets_assembly_fasta 
    String reference_fetch_ncbi_datasets_version = ncbi_datasets_download_genome_accession.ncbi_datasets_version
    String reference_fetch_ncbi_datasets_docker = ncbi_datasets_download_genome_accession.ncbi_datasets_docker
  }

}