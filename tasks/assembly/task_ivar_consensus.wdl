version 1.0

task consensus {
  input {
    File bamfile
    String samplename
    File? reference_genome
    Boolean count_orphans = true
    Int max_depth = "600000"
    Boolean disable_baq = true
    Int min_bq = "0"
    Int min_qual = "20"
    Float? consensus_min_freq 
    Int? consensus_min_depth
    String char_unknown = "N"
    Int disk_size = 100
  }
  command <<<
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
    
    # call consensus
    samtools mpileup \
    ~{true = "--count-orphans" false = "" count_orphans} \
    -d ~{max_depth} \
    ~{true = "--no-BAQ" false = "" disable_baq} \
    -Q ~{min_bq} \
    --reference ${ref_genome} \
    ~{bamfile} | \
    ivar consensus \
    -p ~{samplename}.consensus \
    -q ~{min_qual} \
    -t ~{consensus_min_freq} \
    -m ~{consensus_min_depth} \
    -n ~{char_unknown}

    # clean up fasta header
    echo ">~{samplename}" > ~{samplename}.ivar.consensus.fasta
    grep -v ">" ~{samplename}.consensus.fa >> ~{samplename}.ivar.consensus.fasta
  >>>
  output {
    File consensus_seq = "~{samplename}.ivar.consensus.fasta"
    String ivar_version = read_string("IVAR_VERSION")
    String samtools_version = read_string("SAMTOOLS_VERSION")
    String pipeline_date = read_string("DATE")
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/staphb/ivar:1.3.1-titan"
    memory: "8 GB"
    cpu: 2
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}