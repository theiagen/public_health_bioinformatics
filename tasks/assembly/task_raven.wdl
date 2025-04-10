version 1.0

task raven {
  input {
    File read1
    String samplename

    Int raven_polishing_iterations = 2
    Float? raven_identity = 0
    String? additional_parameters # Any extra Raven-specific parameters

    Int cpu = 4
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/raven:1.8.3"
    Int memory = 16
  }
  command <<<
    # Fail hard
    set -euo pipefail

    # Date and version control
    raven --version | tee VERSION

    # Run Raven 
    raven \
      --polishing-rounds ~{raven_polishing_iterations} \
      ~{"--identity " + raven_identity} \
      ~{"--extra-params " + additional_parameters } \
      --threads ~{cpu} \
      ~{read1} > ~{samplename}.assembly.fasta

    # Check if assembly was successful by checking for output file content
    if [ ! -s ~{samplename}.assembly.fasta ]; then
      echo "Raven assembly failed. No output file generated."
      exit 1
    fi
  >>>
  output {
    File assembly_fasta = "~{samplename}.assembly.fasta"
    String raven_version = read_string("VERSION")
    String raven_docker = "~{docker}"
  }
  runtime {
    docker: "~{docker}"
    cpu: cpu
    memory: "~{memory} GB"
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}