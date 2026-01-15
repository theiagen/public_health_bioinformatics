version 1.0

import "../../../tasks/phylogenetic_inference/utilities/task_referenceseeker.wdl" as referenceseeker_task
import "../../../tasks/task_versioning.wdl" as versioning
import "../../../tasks/utilities/data_import/task_ncbi_datasets.wdl" as ncbi_datasets_task

workflow assembly_fetch {
  input {
    String samplename
    File? assembly_fasta
    String? ncbi_accession
  }
  # if user does not provide reference genome, determine one for the user by running referenceseeker and ncbi datasets on a provided genome to acquire one
  if(! defined(ncbi_accession) && defined(assembly_fasta)) {
    call referenceseeker_task.referenceseeker {
      input:
        assembly_fasta = select_first([assembly_fasta, ""]), # to avoid error with required task input
        samplename = samplename
    }
  }
  call ncbi_datasets_task.ncbi_datasets_download_genome_accession {
    input:
      ncbi_accession = select_first([ncbi_accession, referenceseeker.referenceseeker_top_hit_ncbi_accession])
  }
  call versioning.version_capture {
    input:
  }
  output {
    # version capture
    String assembly_fetch_version = version_capture.phb_version
    String assembly_fetch_analysis_date = version_capture.date
    # referenceseeker outputs
    String? assembly_fetch_referenceseeker_top_hit_ncbi_accession = referenceseeker.referenceseeker_top_hit_ncbi_accession
    String? assembly_fetch_referenceseeker_version = referenceseeker.referenceseeker_version
    File? assembly_fetch_referenceseeker_tsv = referenceseeker.referenceseeker_tsv
    String? assembly_fetch_referenceseeker_docker = referenceseeker.referenceseeker_docker
    String? assembly_fetch_referenceseeker_database = referenceseeker.referenceseeker_database
    # ncbi datasets outputs
    File? assembly_fetch_ncbi_datasets_assembly_data_report_json = ncbi_datasets_download_genome_accession.ncbi_datasets_assembly_data_report_json
    File? assembly_fetch_ncbi_datasets_assembly_fasta = ncbi_datasets_download_genome_accession.ncbi_datasets_assembly_fasta 
    File? assembly_fetch_ncbi_datasets_gff3 = ncbi_datasets_download_genome_accession.ncbi_datasets_gff3
    File? assembly_fetch_ncbi_datasets_gff = ncbi_datasets_download_genome_accession.ncbi_datasets_gbff
    String assembly_fetch_ncbi_datasets_version = ncbi_datasets_download_genome_accession.ncbi_datasets_version
    String assembly_fetch_ncbi_datasets_docker = ncbi_datasets_download_genome_accession.ncbi_datasets_docker
  }

}