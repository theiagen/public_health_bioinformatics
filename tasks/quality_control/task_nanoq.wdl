version 1.0

task nanoq {
  input {
    File reads
    String samplename
    String docker = "quay.io/biocontainers/nanoq:0.9.0--hec16e2b_1"
    Int disk_size = 100
    Int max_read_length = 10000
    Int min_read_length = 500
    Int max_read_qual = 40
    Int min_read_qual = 10
  }
  command <<<
    # capture date and version
    nanoq --version | grep nanoq | tee VERSION

    nanoq -i ~{reads} --min-len ~{min_read_length} --max-len ~{max_read_length} --min-qual ~{min_read_qual} --max-qual ~{max_read_qual} -o ~{samplename}_reads.fastq.gz
  >>>
  output {
    File filtered_reads = "${samplename}_reads.fq.gz"
    String version = read_string("VERSION")
  }
  runtime {
    docker:  "~{docker}"
    memory:  "2 GB"
    cpu:   2
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible:  0
  }
}