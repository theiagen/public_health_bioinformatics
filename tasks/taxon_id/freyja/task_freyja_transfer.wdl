version 1.0

task transfer_files {
  input {
    String gcp_uri
    File updated_barcodes
    File updated_lineages
    File update_log
    Int disk_size = 25
    Int memory = 2
    Int cpu = 1
    String docker = "us-docker.pkg.dev/general-theiagen/cloudsdktool/google-cloud-cli:427.0.0-alpine"
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