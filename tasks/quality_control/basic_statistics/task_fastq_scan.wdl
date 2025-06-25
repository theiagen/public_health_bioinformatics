version 1.0

task fastq_scan_pe {
  input {
    File read1
    File read2
    String read1_name = basename(basename(basename(read1, ".gz"), ".fastq"), ".fq")
    String read2_name = basename(basename(basename(read2, ".gz"), ".fastq"), ".fq")
    Int disk_size = 50
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/fastq-scan:1.0.1--h4ac6f70_3"
    Int memory = 4
    Int cpu = 1
  }
  command <<<
    # exit task in case anything fails in one-liners or variables are unset
    set -euo pipefail

    # capture version
    fastq-scan -v | tee VERSION

    # set cat command based on compression
    if [[ "~{read1}" == *".gz" ]] ; then
      cat_reads="zcat"
    else
      cat_reads="cat"
    fi

    # capture forward read stats
    echo "DEBUG: running fastq-scan on $(basename ~{read1})"
    eval "${cat_reads} ~{read1}" | fastq-scan | tee ~{read1_name}_fastq-scan.json
    # using simple redirect so STDOUT is not confusing
    jq .qc_stats.read_total ~{read1_name}_fastq-scan.json > READ1_SEQS
    echo "DEBUG: number of reads in $(basename ~{read1}): $(cat READ1_SEQS)"
    read1_seqs=$(cat READ1_SEQS)
    echo

    # capture reverse read stats
    echo "DEBUG: running fastq-scan on $(basename ~{read2})"
    eval "${cat_reads} ~{read2}" | fastq-scan | tee ~{read2_name}_fastq-scan.json

    # using simple redirect so STDOUT is not confusing
    jq .qc_stats.read_total ~{read2_name}_fastq-scan.json > READ2_SEQS
    echo "DEBUG: number of reads in $(basename ~{read2}): $(cat READ2_SEQS)"
    read2_seqs=$(cat READ2_SEQS)

    # capture number of read pairs
    if [ "${read1_seqs}" == "${read2_seqs}" ]; then
      read_pairs=${read1_seqs}
    else
      read_pairs="Uneven pairs: R1=${read1_seqs}, R2=${read2_seqs}"
    fi
    
    # use simple redirect so STDOUT is not confusing
    echo "$read_pairs" > READ_PAIRS
    echo "DEBUG: number of read pairs: $(cat READ_PAIRS)"
  >>>
  output {
    File read1_fastq_scan_json = "~{read1_name}_fastq-scan.json"
    File read2_fastq_scan_json = "~{read2_name}_fastq-scan.json"
    Int read1_seq = read_int("READ1_SEQS")
    Int read2_seq = read_int("READ2_SEQS")
    String read_pairs = read_string("READ_PAIRS")
    String version = read_string("VERSION")
    String fastq_scan_docker = docker
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 1
    maxRetries: 3
  }
}

task fastq_scan_se {
  input {
    File read1
    String read1_name = basename(basename(basename(read1, ".gz"), ".fastq"), ".fq")
    Int disk_size = 50
    Int memory = 2
    Int cpu = 1
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/fastq-scan:1.0.1--h4ac6f70_3"
  }
  command <<<
    # exit task in case anything fails in one-liners or variables are unset
    set -euo pipefail
    
    # capture version
    fastq-scan -v | tee VERSION

    # set cat command based on compression
    if [[ "~{read1}" == *".gz" ]] ; then
      cat_reads="zcat"
    else
      cat_reads="cat"
    fi

    # capture forward read stats
    echo "DEBUG: running fastq-scan on $(basename ~{read1})"
    eval "${cat_reads} ~{read1}" | fastq-scan | tee ~{read1_name}_fastq-scan.json
    # using simple redirect so STDOUT is not confusing
    jq .qc_stats.read_total ~{read1_name}_fastq-scan.json > READ1_SEQS
    echo "DEBUG: number of reads in $(basename ~{read1}): $(cat READ1_SEQS)"
  >>>
  output {
    File fastq_scan_json = "~{read1_name}_fastq-scan.json"
    Int read1_seq = read_int("READ1_SEQS")
    String version = read_string("VERSION")
    String fastq_scan_docker = docker
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 1
    maxRetries: 3
  }
}
