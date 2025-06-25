version 1.0

task raven {
  input {
    File read1
    String samplename

    Int raven_polishing_iterations = 2
    Float? raven_identity = 0
    String? raven_opts # Any extra Raven-specific parameters

    Int cpu = 4
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/raven:1.8.3"
    Int memory = 16
  }
  command <<<
    # date and version control
    raven --version | tee VERSION

    # run Raven 
    raven \
      --polishing-rounds ~{raven_polishing_iterations} \
      ~{"--identity " + raven_identity} \
      ~{"--extra-params " + raven_opts} \
      --threads ~{cpu} \
      ~{read1} > ~{samplename}_contigs.fasta

    # check if assembly was successful by checking for output file content
    if [ ! -s ~{samplename}_contigs.fasta ]; then
      echo "DEBUG: Raven assembly failed. No output generated."
      echo "FAIL" > STATUS
    else
      echo "DEBUG: Raven assembly completed successfully."
      echo "PASS" > STATUS
    fi
  >>>
  output {
    File? assembly_fasta = "~{samplename}_contigs.fasta"
    String raven_status = read_string("STATUS")
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