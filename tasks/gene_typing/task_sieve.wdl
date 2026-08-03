version 1.0

task sieve {
  input {
    String samplename
    String plugin
    File? assembly
    File? read1
    File? read2

    File? database # tar.gz compressed database
    String? engine # blast / kma
    String output_format = "text" # text / json / tsv / csv

    Float? min_identity
    Float? min_coverage

    String? parameters # comma-delimited list of parameters to populate

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/sieve:v0.1.0"
    Int disk_size = 50
    Int cpu = 1
    Int memory = 4
  }
  command <<<
  


  >>>
  output {
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 1 # does not take long (usually <3 min) to run stxtyper on 1 genome, preemptible is fine
    maxRetries: 3
  }
}
