version 1.0

task chroquetas {
  input {
    File input_fasta
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/chroquetas:1.0.0"
    Int cpu = 2
    Int disk_size = 20
    Int memory = 8
  }
  command <<< 
  >>>
  output {
  }

  runtime {
    memory: "~{memory} GB"
    cpu: cpu
    docker: docker
    disks: "local-disk " + disk_size + " SSD"
    maxRetries: 1
  }
}
