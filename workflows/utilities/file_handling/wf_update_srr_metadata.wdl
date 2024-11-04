version 1.0

import "../../../tasks/utilities/data_handling/task_fetch_srr_metadata.wdl" as srr_task

workflow wf_retrieve_srr {
  meta {
    description: "This workflow retrieves the Sequence Read Archive (SRA) accession (SRR) associated with a given BioSample accession. It uses the fastq-dl tool to fetch metadata from SRA and outputs the SRR accession that can be used for downstream analysis."
  }
  input {
    String biosample_accession  
    String retrieve_srr_docker = "us-docker.pkg.dev/general-theiagen/biocontainers/fastq-dl:2.0.4--pyhdfd78af_0"
  }

  call srr_task.get_srr_from_biosamples {
    input:
      biosample_accession = biosample_accession,  
      docker = retrieve_srr_docker
  }

  output {
    String biosample_srr_accession = get_srr_from_biosamples.srr_accession 
  }
}
