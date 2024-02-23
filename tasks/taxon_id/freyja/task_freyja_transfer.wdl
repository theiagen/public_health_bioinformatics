version 1.0

task transfer_files {
  input {
    String gcp_uri
    File updated_barcodes
    File updated_lineages
    File update_log
    Int disk_size = 100
    Int memory = 4
    Int cpu = 2
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1"
  }
  command <<<
  # transfer_files to specified gcp_uri
  date_tag=$(date +"%Y-%m-%d")
  
  gsutil -m cp ~{updated_barcodes} ~{updated_lineages} ~{update_log} ~{gcp_uri}/${date_tag}
  
  >>>
  runtime {
    memory: memory + " GB"
    cpu: cpu
    docker: docker
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
  }
  output {    
  }
}