version 1.0

task zip_files {
  input {
    Array[File] files_to_zip
    String zipped_file_name
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
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}