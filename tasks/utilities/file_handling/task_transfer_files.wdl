version 1.0

task transfer_files {
  input {
    Array[String] files_to_transfer
    String target_bucket
    Int cpus = 4
    Int mem_size_gb = 8
    String docker_image = "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1"
  }
  command <<<
    file_path_array="~{sep=' ' files_to_transfer}"

    gsutil -m cp -n ${file_path_array[@]} ~{target_bucket}
    
    echo "transferred_files" > transferred_files.tsv
    gsutil ls ~{target_bucket} >> transferred_files.tsv        
   >>>
  output {
    File transferred_files = "transferred_files.tsv"
  }
  runtime {
      docker: "~{docker_image}"
      memory: "~{mem_size_gb} GB"
      cpu: cpus
      disks: "local-disk 100 SSD"
      preemptible: 0
  }
}