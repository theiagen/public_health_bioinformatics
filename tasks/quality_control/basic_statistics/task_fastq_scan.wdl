version 1.0

task fastq_scan_pe {
  input {
    File read1
    File read2
    Int disk_size = 50
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/fastq-scan:1.0.1--h4ac6f70_3"
    Int memory = 4
    Int cpu = 1
  }
  String read1_name = basename(basename(basename(read1, ".gz"), ".fastq"), ".fq")
  String read2_name = basename(basename(basename(read2, ".gz"), ".fastq"), ".fq")
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

    jq -r .qc_stats.read_mean ~{read1_name}_fastq-scan.json > READ1_MEAN_LENGTH
    echo "DEBUG: mean read length in $(basename ~{read1}): $(cat READ1_MEAN_LENGTH)"

    jq -r .qc_stats.qual_mean ~{read1_name}_fastq-scan.json > READ1_MEAN_QUALITY
    echo "DEBUG: mean read quality in $(basename ~{read1}): $(cat READ1_MEAN_QUALITY)"

    # capture reverse read stats
    echo "DEBUG: running fastq-scan on $(basename ~{read2})"
    eval "${cat_reads} ~{read2}" | fastq-scan | tee ~{read2_name}_fastq-scan.json

    # using simple redirect so STDOUT is not confusing
    jq .qc_stats.read_total ~{read2_name}_fastq-scan.json > READ2_SEQS
    echo "DEBUG: number of reads in $(basename ~{read2}): $(cat READ2_SEQS)"
    read2_seqs=$(cat READ2_SEQS)

    jq -r .qc_stats.read_mean ~{read2_name}_fastq-scan.json > READ2_MEAN_LENGTH
    echo "DEBUG: mean read length in $(basename ~{read2}): $(cat READ2_MEAN_LENGTH)"

    jq -r .qc_stats.qual_mean ~{read2_name}_fastq-scan.json > READ2_MEAN_QUALITY
    echo "DEBUG: mean read quality in $(basename ~{read2}): $(cat READ2_MEAN_QUALITY)"

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
    Float read1_mean_length = read_float("READ1_MEAN_LENGTH")
    Float read2_mean_length = read_float("READ2_MEAN_LENGTH")
    Float read1_mean_quality = read_float("READ1_MEAN_QUALITY")
    Float read2_mean_quality = read_float("READ2_MEAN_QUALITY")
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
    Int disk_size = 50
    Int memory = 2
    Int cpu = 1
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/fastq-scan:1.0.1--h4ac6f70_3"
  }  
  String read1_name = basename(basename(basename(read1, ".gz"), ".fastq"), ".fq")
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

    jq -r .qc_stats.read_mean ~{read1_name}_fastq-scan.json > READ1_MEAN_LENGTH
    echo "DEBUG: mean read length in $(basename ~{read1}): $(cat READ1_MEAN_LENGTH)"

    jq -r .qc_stats.qual_mean ~{read1_name}_fastq-scan.json > READ1_MEAN_QUALITY
    echo "DEBUG: mean read quality in $(basename ~{read1}): $(cat READ1_MEAN_QUALITY)"
  >>>
  output {
    File fastq_scan_json = "~{read1_name}_fastq-scan.json"
    Int read1_seq = read_int("READ1_SEQS")
    Float read1_mean_length = read_float("READ1_MEAN_LENGTH")
    Float read1_mean_quality = read_float("READ1_MEAN_QUALITY")
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
