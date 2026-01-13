version 1.0

task fastq_scan {
  input {
    File read1
    File? read2
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

    # set default values
    echo "0" | tee READ1_SEQS READ2_SEQS
    echo "0.0" | tee READ1_MEAN_LENGTH READ2_MEAN_LENGTH
    echo "0.0" | tee READ1_MEAN_QUALITY READ2_MEAN_QUALITY

    # processes one or both reads as needed
    all_reads="~{read1} ~{read2}"
    for read_file in ${all_reads}; do
      if [[ "$read_file" == "~{read1}" ]]; then
        read_name="~{read1_name}"
        read_label="READ1"
      else
        read_name="~{read2_name}"
        read_label="READ2"
      fi

      # set cat command based on compression
      if [[ "$read_file" == *".gz" ]] ; then
        cat_reads="zcat"
      else
        cat_reads="cat"
      fi

      # capture read stats
      echo "DEBUG: running fastq-scan on $(basename $read_file)"
      eval "${cat_reads} $read_file" | fastq-scan | tee ${read_name}_fastq-scan.json

      # fail gracefully if the output json is empty or missing (e.g., from an empty input fastq)
      if [ ! -s "${read_name}_fastq-scan.json" ]; then
        echo "ERROR: fastq-scan output JSON is missing or empty for $(basename $read_file)."
      else
        # using simple redirect so STDOUT is not confusing
        jq .qc_stats.read_total ${read_name}_fastq-scan.json > ${read_label}_SEQS
        echo "DEBUG: number of reads in $(basename $read_file): $(cat ${read_label}_SEQS)"

        jq -r .qc_stats.read_mean ${read_name}_fastq-scan.json > ${read_label}_MEAN_LENGTH
        echo "DEBUG: mean read length in $(basename $read_file): $(cat ${read_label}_MEAN_LENGTH)"

        jq -r .qc_stats.qual_mean ${read_name}_fastq-scan.json > ${read_label}_MEAN_QUALITY
        echo "DEBUG: mean read quality in $(basename $read_file): $(cat ${read_label}_MEAN_QUALITY)"
      fi
    done

    # capture number of read pairs, if paired-end
    if [ -n "~{read2}" ]; then
      if [ "$(cat READ1_SEQS)" == "$(cat READ2_SEQS)" ]; then
        read_pairs=$(cat READ1_SEQS)
      else
        read_pairs="Uneven pairs: R1=$(cat READ1_SEQS), R2=$(cat READ2_SEQS)"
      fi
    else
      read_pairs="N/A - single-end data"
    fi
    
    # use simple redirect so STDOUT is not confusing
    echo "$read_pairs" > READ_PAIRS
    echo "DEBUG: number of read pairs: $(cat READ_PAIRS)"
  >>>
  output {
    File read1_fastq_scan_json = "~{read1_name}_fastq-scan.json"
    File? read2_fastq_scan_json = "~{read2_name}_fastq-scan.json"
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
