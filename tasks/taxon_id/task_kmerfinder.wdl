version 1.0

task kmerfinder_bacteria {
  input {
    File assembly
    File kmerfinder_db = "gs://theiagen-public-files-rp/terra/theiaprok-files/kmerfinder_bacteria.tar.gz"
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/kmerfinder:3.0.2--hdfd78af_0"
    Int memory = 32
    Int cpu = 4
    String kmerfinder_args = ""
  }
  command <<<
    # Decompress the kmerfinder database
    mkdir db
    tar -C ./db/ -xzvf ~{kmerfinder_db}  

    # Run kmerfinder
    kmerfinder.py \
        -db ./db/bacteria/bacteria.ATG \
        -tax ./db/bacteria/bacteria.tax \
        -i ~{assembly} \
        -o ~{samplename} \
        ~{kmerfinder_args} \
        -x

  >>>
  output {
    String kmerfinder_docker = docker
  }
  runtime {
    docker: docker
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk 100 SSD"
    preemptible: 0
  }
}