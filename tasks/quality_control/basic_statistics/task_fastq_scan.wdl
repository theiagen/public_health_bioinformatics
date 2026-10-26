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
    # Not adding these to JSON because it is currently beyond scope
    echo "DEBUG: running seqkit stats on $(basename -a ~{read1} ~{read2} | paste -sd ' ')"
    seqkit stats -a -T -j ~{cpu} ~{read1} ~{read2} | tee seqkit_stats.tsv
    # column 18 is GC(%) in seqkit v2.14.0 `stats -a -T` output; revisit if the docker's seqkit version changes
    awk -F'\t' 'NR==2 {print $18}' seqkit_stats.tsv > READ1_GC_PERCENT
    echo "DEBUG: GC content in $(basename ~{read1}): $(cat READ1_GC_PERCENT)"
    if [[ -n "~{read2}" ]]; then
      awk -F'\t' 'NR==3 {print $18}' seqkit_stats.tsv > READ2_GC_PERCENT
      echo "DEBUG: GC content in $(basename ~{read2}): $(cat READ2_GC_PERCENT)"
    fi

    # capture mean quality across all reads, weighted by the number of reads in each file
    if [[ -n "~{read2}" ]]; then
      read2_mean_quality=$(cat READ2_MEAN_QUALITY)
    else
      read2_seqs=0
      read2_mean_quality=0
    fi
    awk -v n1="${read1_seqs}" -v v1="$(cat READ1_MEAN_QUALITY)" -v n2="${read2_seqs}" -v v2="${read2_mean_quality}" \
      'BEGIN {n = n1 + n2; if (n > 0) print (n1 * v1 + n2 * v2) / n; else print 0}' > COMBINED_MEAN_QUALITY
    echo "DEBUG: mean read quality across all reads: $(cat COMBINED_MEAN_QUALITY)"

    # capture GC content across all reads, weighted by the number of bases in each file
    # column 5 is sum_len and column 18 is GC(%) in seqkit v2.14.0
    awk -F'\t' 'NR>1 {bases += $5; gc_bases += $5 * $18}
      END {if (bases > 0) print gc_bases / bases; else print 0}' seqkit_stats.tsv > COMBINED_GC_PERCENT
    echo "DEBUG: GC content across all reads: $(cat COMBINED_GC_PERCENT)"
  >>>
  output {
    File read1_fastq_scan_json = "~{read1_name}_fastq-scan.json"
    File? read2_fastq_scan_json = "~{read2_name}_fastq-scan.json"
    Int read1_seq = read_int("READ1_SEQS")
    # 0 Values for non-existent read2 because they are not propagated in SE/ONT tasks
    Int read2_seq = if defined(read2) then read_int("READ2_SEQS") else 0
    Float read1_mean_length = read_float("READ1_MEAN_LENGTH")
    Float read2_mean_length = if defined(read2) then read_float("READ2_MEAN_LENGTH") else 0.0
    Float read1_mean_quality = read_float("READ1_MEAN_QUALITY")
    Float read2_mean_quality = if defined(read2) then read_float("READ2_MEAN_QUALITY") else 0.0
    Float read1_gc_percent = read_float("READ1_GC_PERCENT")
    Float read2_gc_percent = if defined(read2) then read_float("READ2_GC_PERCENT") else 0.0
    Float total_mean_quality = read_float("COMBINED_MEAN_QUALITY")
    Float total_gc_percent = read_float("COMBINED_GC_PERCENT")
    String read_pairs = if defined(read2) then read_string("READ_PAIRS") else "NA"
    String fastq_scan_version = read_string("VERSION")
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
