version 1.0

task consensus {
  input {
    File bamfile 
    String samplename
    File? reference_genome
    Boolean count_orphans = true
    Int max_depth = "600000"
    Boolean disable_baq = true
    Boolean all_positions = false
    Int min_bq = "0"
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
    samtools --version | head -n1 | tee SAMTOOLS_VERSION

    # set reference genome
    if [[ ! -z "~{reference_genome}" ]]; then
      echo "User reference identified; ~{reference_genome} will be utilized for alignement"
      ref_genome="~{reference_genome}"
      bwa index "~{reference_genome}"
      # move to primer_schemes dir; bwa fails if reference file not in this location
    else
      ref_genome="/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta"  
    fi

    # call variants
    samtools mpileup \
      ~{true = "-A" false = "" count_orphans} \
      -d ~{max_depth} \
      ~{true = "-B" false = "" disable_baq} \
      -Q ~{min_bq} \
      --reference ${ref_genome} \
      ~{true = "-aa" false = "" all_positions} \
      ~{bamfile} \
      > ~{samplename}.mpileup

    cat ~{samplename}.mpileup | \
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
    File sample_mpileup = "~{samplename}.mpileup"
    String ivar_version = read_string("IVAR_VERSION")
    String pipeline_date = read_string("DATE")
    String samtools_version = read_string("SAMTOOLS_VERSION")
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
