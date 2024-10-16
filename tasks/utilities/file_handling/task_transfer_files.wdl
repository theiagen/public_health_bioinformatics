version 1.0

task transfer_files {
  input {
    Array[String] files_to_transfer
    String target_bucket
    Int cpu = 4
    Int memory = 8
    String docker_image = "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1"
    Int disk_size = 100
  }
  meta {
    # added so that call caching is always turned off
    volatile: true
  }
  command <<<
    set -euo pipefail
    
    file_path_array="~{sep=' ' files_to_transfer}"

    gsutil -m cp -n ${file_path_array[@]} ~{target_bucket}
    
    echo "transferred_files" > transferred_files.tsv
    gsutil -m ls ~{target_bucket} >> transferred_files.tsv        
   >>>
  output {
    File transferred_files = "transferred_files.tsv"
  }
  runtime {
    docker: "~{docker_image}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}