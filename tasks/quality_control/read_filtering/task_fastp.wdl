version 1.0

task fastp_pe {
  input {
    File read1
    File read2
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/fastp:0.23.2"
    Int disk_size = 100
    Int fastp_window_size = 4 # set to mirror v0.23.2 default 
    Int fastp_quality_trim_score = 20 # set to mirror v0.23.2 default
    Int fastp_min_length = 15 # set to mirror v0.23.2 default
    # -g enables polyg trimming with default value of 10
    String? fastp_args
    Int cpu = 4
    Int memory = 8
  }
  command <<<
    # fail hard
    set -euo pipefail

    # date 
    date | tee DATE

    # version
    fastp -v 2>&1 | sed -E 's/^fastp //g' | tee VERSION

    # trim reads
    fastp \
      --in1 ~{read1} --in2 ~{read2} \
      --out1 ~{samplename}_1P.fastq.gz --out2 ~{samplename}_2P.fastq.gz \
      --unpaired1 ~{samplename}_1U.fastq.gz --unpaired2 ~{samplename}_2U.fastq.gz \
      ~{if defined(fastp_window_size) || defined(fastp_quality_trim_score) || defined(fastp_min_length) then "--cut_right" else ""} \
      ~{if defined(fastp_window_size) then "--cut_right_window_size ~{fastp_window_size}" else ""} \
      ~{if defined(fastp_quality_trim_score) then "--cut_right_mean_quality ~{fastp_quality_trim_score}" else ""} \
      ~{if defined(fastp_min_length) then "--length_required ~{fastp_min_length}" else ""} \
      --thread ~{cpu} \
      ~{if defined(fastp_args) then "~{fastp_args}" else ""} \
      --html ~{samplename}_fastp.html --json ~{samplename}_fastp.json
  >>>
  output {
    File read1_trimmed = "~{samplename}_1P.fastq.gz"
    File read2_trimmed = "~{samplename}_2P.fastq.gz"
    File read1_trimmed_unpaired = "~{samplename}_1U.fastq.gz"
    File read2_trimmed_unpaired = "~{samplename}_2U.fastq.gz"
    File fastp_stats_html = "~{samplename}_fastp.html"
    File fastp_stats_json = "~{samplename}_fastp.json"
    String fastp_version = read_string("VERSION")
    String fastp_docker = "~{docker}"
    String pipeline_date = read_string("DATE")
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}

task fastp_se {
  input {
    File read1
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/fastp:0.23.2"
    Int disk_size = 100
    Int fastp_window_size = 20
    Int fastp_quality_trim_score = 30
    Int fastp_min_length = 50
    # -g enables polyg trimming with default value of 10
    # --detect_adapter_for_pe argument was removed 
    String fastp_args = "-g -5 20 -3 20"
    Int cpu = 4
    Int memory = 8
  }
  command <<<
    # date 
    date | tee DATE

    fastp \
    --in1 ~{read1} \
    --out1 ~{samplename}_1P.fastq.gz \
    --cut_right --cut_right_window_size ~{fastp_window_size} --cut_right_mean_quality ~{fastp_quality_trim_score} \
    --length_required ~{fastp_min_length} \
    --thread ~{cpu} \
    ~{fastp_args} \
    --html ~{samplename}_fastp.html --json ~{samplename}_fastp.json
  >>>
  output {
    File read1_trimmed = "~{samplename}_1P.fastq.gz"
    File fastp_stats = "~{samplename}_fastp.html"
    String version = "~{docker}"
    String pipeline_date = read_string("DATE")
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}