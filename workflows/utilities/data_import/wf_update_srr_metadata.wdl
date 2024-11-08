version 1.0

import "../../../tasks/utilities/data_handling/task_fetch_srr_metadata.wdl" as srr_task

workflow wf_retrieve_srr {
  meta {
    description: "This workflow retrieves the Sequence Read Archive (SRA) accession (SRR) associated with a given sample accession. It uses the fastq-dl tool to fetch metadata from SRA and outputs the SRR accession that can be used for downstream analysis."
  }
  input {
    String sample_accession  
  }

  call srr_task.fetch_srr_metadata {
    input:
      sample_accession = sample_accession
  }

  output {
    String srr_accession = fetch_srr_metadata.srr_accession 
  }
}