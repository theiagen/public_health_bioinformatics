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
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/augur:31.5.0"
  }
  command <<<
    set -euo pipefail
    
    # capture version information
    augur version > VERSION
    echo
    echo "mafft version:"
    mafft --version 2>&1 | tee MAFFT_VERSION

    # run augur align
    augur align \
      --sequences ~{assembly_fasta} \
      --reference-sequence ~{reference_fasta} \
      --nthreads ~{cpu} \
      ~{true="--remove-reference" false="" remove_reference} \
      ~{true="--fill-gaps" false="" fill_gaps}
  >>>
  output {
    File aligned_fasta = "alignment.fasta"
    String augur_version = read_string("VERSION")
    String mafft_version = read_string("MAFFT_VERSION")
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