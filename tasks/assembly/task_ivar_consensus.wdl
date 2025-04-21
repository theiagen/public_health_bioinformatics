version 1.0

task consensus {
  input {
    File mpileup
    String samplename
    Int min_qual = "20"
    Float? consensus_min_freq 
    Int? consensus_min_depth
    String char_unknown = "N"
    Boolean skip_N = false
    Int disk_size = 100
    Int memory = 8
    Int cpu = 2
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/ivar:1.3.1-titan"
  }
  command <<<
    #set -euo pipefail to avoid silent failure
    set -euo pipefail
    # date and version control
    date | tee DATE
    ivar version | head -n1 | tee IVAR_VERSION

    cat ~{mpileup} | \
    ivar consensus \
      -p ~{samplename}.consensus \
      -q ~{min_qual} \
      -t ~{consensus_min_freq} \
      -m ~{consensus_min_depth} \
      -n ~{char_unknown} \
      ~{true = "-k" false = "" skip_N} 

    # clean up fasta header
    echo ">~{samplename}" > ~{samplename}.ivar.consensus.fasta
    grep -v ">" ~{samplename}.consensus.fa >> ~{samplename}.ivar.consensus.fasta
  >>>
  output {
    File consensus_seq = "~{samplename}.ivar.consensus.fasta"
    String ivar_version = read_string("IVAR_VERSION")
    String pipeline_date = read_string("DATE")
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
