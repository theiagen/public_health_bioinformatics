version 1.0

import "../../../tasks/utilities/data_handling/task_fetch_srr_accession.wdl" as srr_task

workflow wf_retrieve_srr {
  meta {
    description: "This workflow retrieves the Sequence Read Archive (SRA) accession (SRR) associated with a given sample accession. It uses the fastq-dl tool to fetch metadata from SRA and outputs the SRR accession."
  }
  input {
    String sample_accession  
  }
  call versioning_task.version_capture {
    input:
  }
  call srr_task.fetch_srr_accession {
    input:
      sample_accession = sample_accession
  }
  output {
    Array[String] srr_accession = fetch_srr_accession.srr_accession
    
    # Version Captures
    String phb_version = version_capture.phb_version
    String fetch_srr_date = version_capture.date
  }
}