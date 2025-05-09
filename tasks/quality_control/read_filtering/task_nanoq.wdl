version 1.0

task nanoq {
  input {
    File read1 # intended for ONT data only
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/nanoq:0.9.0--hec16e2b_1"
    Int disk_size = 100
    Int max_read_length = 100000
    Int min_read_length = 500
    Int max_read_qual = 100 # set arbitrarily high to capture all reads
    Int min_read_qual = 10
    Int memory = 2
    Int cpu = 2
  }
  command <<<
    # capture date and version
    nanoq --version | grep nanoq | tee VERSION

    nanoq -i ~{read1} --min-len ~{min_read_length} --max-len ~{max_read_length} --min-qual ~{min_read_qual} --max-qual ~{max_read_qual} -o ~{samplename}_read1.fastq.gz
  >>>
  output {
    File filtered_read1 = "${samplename}_read1.fastq.gz"
    String version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible:  0
  }
}