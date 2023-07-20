version 1.0

task augur_align {
  input {
    File assembly_fasta
    File reference_fasta
    Boolean fill_gaps = false
    Int cpus = 64
    Int mem_size = 32
    Int disk_size = 750
  }
  command <<<
    # capture version information
    augur version > VERSION

    # run augur align
    augur align \
      --sequences ~{assembly_fasta} \
      --nthreads ~{cpus} \
      --reference-sequence ~{reference_fasta} \
      ~{true="--fill-gaps" false="" fill_gaps}
  >>>
  output {
    File aligned_fasta = "alignment.fasta"
    String augur_version = read_string("VERSION")
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/biocontainers/augur:22.0.2--pyhdfd78af_0"
    memory: mem_size + " GB"
    cpu :   cpus
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    preemptible: 0
    dx_instance_type: "mem3_ssd1_v2_x36"
    maxRetries: 3
  }
}