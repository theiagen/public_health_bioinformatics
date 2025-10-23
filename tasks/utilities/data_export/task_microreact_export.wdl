version 1.0

task microreact_export {
  meta {
    description: "This task generates a Microreact input file from supplied metadata and tree files, submitting to Microreact when access token is provided."
  }
  input {
    File metadata
    Array[String] metadata_columns
    Array[String] tree_files
    String docker 
    Int disk_size = 10
    Int memory = 4
    Int cpu = 2
  }
  command <<<
    # set -euo pipefail to avoid silent failure
    set -euo pipefail

    
  >>>
  output {
    File microreact_json 
    File? microreact_api_response 
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}