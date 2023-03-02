version 1.0

task test {
  input {
    File reference_genome
    String docker_image = "quay.io/staphb/lyveset:1.1.4f"
    Int memory = 16
    Int cpu = 4
    Int disk_size = 100
  }
  command <<<
    date | tee DATE
    echo "tree not created" > tree_not.txt
  >>>
  output {
    String lyveset_docker_image = docker_image
    File ref_fasta = "~{reference_genome}"
    File lyveset_log = stdout()
  }
  runtime {
    docker: docker_image
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk ~{disk_size} SSD"
    preemptible: 0
    maxRetries: 0
  }
}