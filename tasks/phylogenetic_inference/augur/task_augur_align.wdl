version 1.0

task augur_align {
  input {
    File assembly_fasta
    File reference_fasta
    Boolean remove_reference
    Boolean fill_gaps = false
    Int cpu = 64
    Int memory = 128
    Int disk_size = 750
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/augur:22.0.2--pyhdfd78af_0"
  }
  command <<<
    # capture version information
    augur version > VERSION

    # run augur align
    augur align \
      --sequences ~{assembly_fasta} \
      --nthreads cpu \
      --reference-sequence ~{reference_fasta} \
      ~{true="--remove-reference" false="" remove_reference} \
      ~{true="--fill-gaps" false="" fill_gaps}
  >>>
  output {
    File aligned_fasta = "alignment.fasta"
    String augur_version = read_string("VERSION")
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    preemptible: 0
    dx_instance_type: "mem3_ssd1_v2_x36"
    maxRetries: 3
  }
}