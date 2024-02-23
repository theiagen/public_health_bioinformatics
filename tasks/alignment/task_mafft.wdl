version 1.0

task mafft {
  input {
    Array[File] genomes
    Int cpu = 16
    Int memory = 32
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/mafft:7.450"
  }
  command <<<
    # date and version control
    date | tee DATE
    mafft_vers=$(mafft --version)
    echo Mafft $(mafft_vers) | tee VERSION

    cat ~{sep=" " genomes} | sed 's/Consensus_//;s/.consensus_threshold.*//' > assemblies.fasta
    mafft --thread -~{cpu} assemblies.fasta > msa.fasta
  >>>
  output {
    String date = read_string("DATE")
    String version = read_string("VERSION")
    File msa = "msa.fasta"
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}