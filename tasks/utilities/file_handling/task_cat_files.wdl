version 1.0

task cat_files {
  input {
    Array[File] files_to_cat
    String concatenated_file_name
    String docker_image = "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1"
    Boolean skip_extra_headers = false
    Int memory = 8
    Int cpu = 2
    Int disk_size = 100
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
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}

task cat_files_fasta {
  input {
    Array[File] files_to_cat
    Array[String] headers
    String concatenated_file_name
    String docker_image = "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1"
    Int memory = 8
    Int cpu = 2
    Int disk_size = 100
  }
  meta {
    # added so that call caching is always turned off
    volatile: true
  }
  command <<<
    file_array=(~{sep=' ' files_to_cat})
    headers_array=(~{sep=' ' headers})
    touch ~{concatenated_file_name}

    # cat files one by one and store them in the concatenated_files file
    for index in ${!file_array[@]}; do
      file=${file_array[$index]}
      header=${headers_array[$index]}
      # replace the original header with the new header using sed before concatenating
      awk 1 ${file} | sed "s/^>.*/>${header}/" >> ~{concatenated_file_name}
    done
  >>>
  output {
    File concatenated_files = "~{concatenated_file_name}"
  }
  runtime {
    docker: "~{docker_image}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}