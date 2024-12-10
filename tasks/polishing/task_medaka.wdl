version 1.0

task medaka_consensus {
  input {
    File unpolished_fasta
    String samplename
    File read1
    String medaka_model = "r1041_e82_400bps_sup_v5.0.0"
    Int cpu = 4
    Int memory = 16
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/medaka:2.0.1"
  }
  command <<< 
    set -euo pipefail

    medaka --version | tee MEDAKA_VERSION

    # Perform a single round of Medaka polishing
    medaka_consensus \
      -i ~{read1} \
      -d ~{unpolished_fasta} \
      -o . \
      -m ~{medaka_model} \
      -t ~{cpu}

    # Rename final polished output with sample name
    mv consensus.fasta ~{samplename}.polished.fasta
  >>>
  output {
    File medaka_fasta = "~{samplename}.polished.fasta"
    String medaka_version = read_string("MEDAKA_VERSION")
  }
  runtime {
    docker: "~{docker}"
    cpu: cpu
    memory: "~{memory} GB"
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    maxRetries: 1
    preemptible: 0
  }
}
