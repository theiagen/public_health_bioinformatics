version 1.0

task transfer_files {
  input {
    String gcp_uri
    File updated_barcodes
    File updated_lineages
    File update_log
    Int disk_size = 100
  }
  command <<<
  # transfer_files to specified gcp_uri
  date_tag=$(date +"%Y-%m-%d")
  
  gsutil -m cp ~{updated_barcodes} ~{updated_lineages} ~{update_log} ~{gcp_uri}/${date_tag}
  
  >>>
  runtime {
    memory: "4 GB"
    cpu: 2
    docker: "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1"
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
  }
  output {    
  }
}