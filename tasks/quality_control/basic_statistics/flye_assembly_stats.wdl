version 1.0

task assembly_stats {
  input {
    File assembly_fasta
    String samplename
    Int cpu = 1
    Int memory = 4
    Int disk_size = 50
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/seqkit:2.8.2"
  }
  command <<< 
    set -euo pipefail

    seqkit stats --all --threads ~{cpu} ~{assembly_fasta} > ~{samplename}_assembly_stats.txt
  >>>
  output {
    File assembly_stats_file = "~{samplename}_assembly_stats.txt"
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}
