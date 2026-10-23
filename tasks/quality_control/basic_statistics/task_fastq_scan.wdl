version 1.0

task fastq_scan {
  input {
    File read1
    File? read2
    Int disk_size = 50
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/fastq-scan:1.0.1-seqkit"
    Int memory = 4
    Int cpu = 1
  }
  String read1_name = basename(basename(basename(read1, ".gz"), ".fastq"), ".fq")
  String read2_name = if defined(read2) then basename(basename(basename(select_first([read2]), ".gz"), ".fastq"), ".fq") else ""
  command <<<
    # exit task in case anything fails in one-liners or variables are unset
    set -euo pipefail

    # capture version
    fastq-scan -v | tee VERSION
    seqkit version | tee SEQKIT_VERSION

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

    if [[ -n "~{read2}" ]]; then
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
    fi

    # capture GC content; seqkit reports one row per file in input order
    echo "DEBUG: running seqkit stats on $(basename -a ~{read1} ~{read2} | paste -sd ' ')"
    seqkit stats -a -T -j ~{cpu} ~{read1} ~{read2} | tee seqkit_stats.tsv
    # locate the GC column by header name since its position differs between seqkit versions
    awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i=="GC(%)") col=i} NR==2 {print $col}' seqkit_stats.tsv > READ1_GC_PERCENT
    echo "DEBUG: GC content in $(basename ~{read1}): $(cat READ1_GC_PERCENT)"
    if [[ -n "~{read2}" ]]; then
      awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i=="GC(%)") col=i} NR==3 {print $col}' seqkit_stats.tsv > READ2_GC_PERCENT
      echo "DEBUG: GC content in $(basename ~{read2}): $(cat READ2_GC_PERCENT)"
    fi
  >>>
  output {
    File read1_fastq_scan_json = "~{read1_name}_fastq-scan.json"
    File? read2_fastq_scan_json = "~{read2_name}_fastq-scan.json"
    Int read1_seq = read_int("READ1_SEQS")
    Int? read2_seq = if defined(read2) then read_int("READ2_SEQS") else 0
    Float read1_mean_length = read_float("READ1_MEAN_LENGTH")
    Float? read2_mean_length = if defined(read2) then read_float("READ2_MEAN_LENGTH") else 0.0
    Float read1_mean_quality = read_float("READ1_MEAN_QUALITY")
    Float? read2_mean_quality = if defined(read2) then read_float("READ2_MEAN_QUALITY") else 0.0
    Float read1_gc_percent = read_float("READ1_GC_PERCENT")
    Float? read2_gc_percent = if defined(read2) then read_float("READ2_GC_PERCENT") else 0.0
    String? read_pairs = if defined(read2) then read_string("READ_PAIRS") else "NA"
    String version = read_string("VERSION")
    String seqkit_version = read_string("SEQKIT_VERSION")
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
