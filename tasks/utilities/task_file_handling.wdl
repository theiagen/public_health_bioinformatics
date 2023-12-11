version 1.0

task cat_files {
  input {
    Array[File] files_to_cat
    String concatenated_file_name
    String docker_image = "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1"
    Boolean skip_extra_headers = false
  }
  meta {
    # added so that call caching is always turned off
    volatile: true
  }
  command <<<
    file_array=(~{sep=' ' files_to_cat})
    touch ~{concatenated_file_name}

    # cat files one by one and store them in the concatenated_files file
    for index in ${!file_array[@]}; do
      file=${file_array[$index]}
      if ! ~{skip_extra_headers} ; then # act as if the first line of all files is not a header
        cat ${file} >> ~{concatenated_file_name}
      else # you want to skip the first line of all files except the first one
        if [ $index == 0 ]; then # if its the first file, cat the entire thing
          cat ${file} >> ~{concatenated_file_name}
        else # otherwise, skip the first line
          tail -n +2 ${file} >> ~{concatenated_file_name}
        fi
      fi
    done
  >>>
  output {
    File concatenated_files = "~{concatenated_file_name}"
  }
  runtime {
    docker: "~{docker_image}"
    memory: "8 GB"
    cpu: 2
    disks: "local-disk 100 SSD"
    preemptible: 0
  }
}

task zip_files {
  input {
    Array[File] files_to_zip
    String zipped_file_name
    String docker_image = "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1"
  }
  meta {
    # added so that call caching is always turned off
    volatile: true
  }
  command <<<
    file_array=(~{sep=' ' files_to_zip})
    mkdir ~{zipped_file_name}

    # move files oto a single directory before zipping
    for index in ${!file_array[@]}; do
      file=${file_array[$index]}
      mv ${file} ~{zipped_file_name}
    done
    
    zip -r ~{zipped_file_name}.zip ~{zipped_file_name}
   >>>
  output {
    File zipped_files = "~{zipped_file_name}.zip"
  }
  runtime {
      docker: "~{docker_image}"
      memory: "8 GB"
      cpu: 2
      disks: "local-disk 100 SSD"
      preemptible: 0
  }
}

task transfer_files {
  input {
    Array[String] files_to_transfer
    String target_bucket
    Int cpus = 4
    Int mem_size_gb = 8
    String docker_image = "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1"
  }
  meta {
    # added so that call caching is always turned off
    volatile: true
  }
  command <<<
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
      memory: "~{mem_size_gb} GB"
      cpu: cpus
      disks: "local-disk 100 SSD"
      preemptible: 0
  }
}