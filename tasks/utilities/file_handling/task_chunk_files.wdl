version 1.0

task chunk_files {
  input {
    File file_list

    Int start_line
    Int end_line

    Int cpu = 1
    Int disk_size = 25
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1"
    Int memory = 2
  }
  command <<<
    sed -n "~{start_line},~{end_line}p" ~{file_list} > chunk.txt
  >>>
  output {
    File chunk = "chunk.txt"
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 1
  }
}
