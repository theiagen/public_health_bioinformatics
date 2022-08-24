version 1.0

task fastq_scan_pe {
  input {
    File read1
    File read2
    String read1_name = basename(basename(basename(read1, ".gz"), ".fastq"), ".fq")
    String read2_name = basename(basename(basename(read2, ".gz"), ".fastq"), ".fq")
  }
  command <<<
    # capture date and version
    date | tee DATE
    fastq-scan -v | tee VERSION

    # set cat command based on compression
    if [[ "~{read1}" == *".gz" ]] ; then
      cat_reads="zcat"
    else
      cat_reads="cat"
    fi

    # capture forward read stats
    eval "${cat_reads} ~{read1}" | fastq-scan | tee ~{read1_name}_fastq-scan.json >(jq .qc_stats.read_total > READ1_SEQS)
    read1_seqs=$(cat READ1_SEQS)
    eval "${cat_reads} ~{read2}" | fastq-scan | tee ~{read2_name}_fastq-scan.json >(jq .qc_stats.read_total > READ2_SEQS)
    read2_seqs=$(cat READ2_SEQS)

    # capture number of read pairs
    if [ "${read1_seqs}" == "${read2_seqs}" ]; then
      read_pairs=${read1_seqs}
    else
      read_pairs="Uneven pairs: R1=${read1_seqs}, R2=${read2_seqs}"
    fi

    echo $read_pairs | tee READ_PAIRS
  >>>
  output {
    File read1_fastq_scan_report = "~{read1_name}_fastq-scan.json"
    File read2_fastq_scan_report = "~{read2_name}_fastq-scan.json"
    Int read1_seq = read_string("READ1_SEQS")
    Int read2_seq = read_string("READ2_SEQS")
    String read_pairs = read_string("READ_PAIRS")
    String version = read_string("VERSION")
    String pipeline_date = read_string("DATE")
  }
  runtime {
    docker: "quay.io/biocontainers/fastq-scan:0.4.4--h7d875b9_1"
    memory: "2 GB"
    cpu: 2
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 3
  }
}

task fastq_scan_se {
  input {
    File read1
    String read1_name = basename(basename(basename(read1, ".gz"), ".fastq"), ".fq")
  }
  command <<<
    # capture date and version
    date | tee DATE
    fastq-scan -v | tee VERSION

    # set cat command based on compression
    if [[ "~{read1}" == *".gz" ]] ; then
      cat_reads="zcat"
    else
      cat_reads="cat"
    fi

    # capture forward read stats
    eval "${cat_reads} ~{read1}" | fastq-scan | tee ~{read1_name}_fastq-scan.json >(jq .qc_stats.read_total > READ1_SEQS)
  >>>
  output {
    File fastq_scan_report = "~{read1_name}_fastq-scan.json"
    Int read1_seq = read_string("READ1_SEQS")
    String version = read_string("VERSION")
    String pipeline_date = read_string("DATE")
  }
  runtime {
    docker: "quay.io/biocontainers/fastq-scan:0.4.4--h7d875b9_1"
    memory: "2 GB"
    cpu: 2
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 3
  }
}