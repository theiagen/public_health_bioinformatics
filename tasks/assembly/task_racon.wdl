version 1.0

task racon {
  input {
    File unpolished_fasta
    File reads
    Int cpu = 4
    Int memory = 16
    Int disk_size = 100
    String samplename
    String docker = "staphb/racon"
  }
  command <<< 
    set -euo pipefail
    minimap2 --version | tee MINIMAP2_VERSION
    racon --version | tee -a RACON_VERSION
    # Generate alignments with Minimap2
    minimap2 -t ~{cpu} ~{unpolished_fasta} ~{reads} > ~{samplename}.paf

    # Run Racon for polishing
    racon \
      -t ~{cpu} \
      ~{reads} \
      ~{samplename}.paf \
      ~{unpolished_fasta} \
      > ~{samplename}.polished.fasta
  >>>
  output {
    File polished_fasta = "~{samplename}.polished.fasta"
    File alignments = "~{samplename}.paf"
    String racon_version = read_string("RACON_VERSION")
    String minimap2_version = read_string("MINIMAP2_VERSION")
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
