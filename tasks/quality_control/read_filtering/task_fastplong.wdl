version 1.0

task fastplong {
  input {
    # required inputs
    File read1
    String samplename

    # quality trimming options
    Int fastplong_window_size = 4 # set to mirror v0.4.1 default 
    Int fastplong_quality_trim_score = 20 # set to mirror v0.4.1 default
    Int fastplong_min_length = 15 # set to mirror v0.4.1 default
    Boolean cut_front = false # 5' to 3'
    Boolean cut_tail = false # 3' to 5'

    # adapter trimming options
    Boolean fastplong_trim_adapters = true
    File? fastplong_adapter_fasta
    String? fastplong_start_adapter
    String? fastplong_end_adapter

    # other options
    String? fastplong_args

    # runtime options
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/fastplong:0.4.1"
    Int disk_size = 100
    Int cpu = 4
    Int memory = 8
  }
  command <<<
    # fail hard
    set -euo pipefail

    # version
    fastplong -v 2>&1 | sed -E 's/^fastplong //g' | tee VERSION

    # trim reads
    fastplong \
      -i ~{read1} \
      -o ~{samplename}_trim.fastq.gz \
      ~{if defined(fastplong_window_size) then "--cut_window_size ~{fastplong_window_size}" else ""} \
      ~{if defined(fastplong_quality_trim_score) then "--cut_mean_quality ~{fastplong_quality_trim_score}" else ""} \
      ~{if cut_front then "--cut_front" else ""} \
      ~{if cut_tail then "--cut_tail" else ""} \
      ~{if defined(fastplong_min_length) then "--length_required ~{fastplong_min_length}" else ""} \
      ~{if ! fastplong_trim_adapters then "--disable_adapter_trimming" else ""} \
      ~{if defined(fastplong_adapter_fasta) then "--adapter_fasta ~{fastplong_adapter_fasta}" else ""} \
      ~{if defined(fastplong_start_adapter) then "--start_adapter ~{fastplong_start_adapter}" else ""} \
      ~{if defined(fastplong_end_adapter) then "--end_adapter ~{fastplong_end_adapter}" else ""} \
      --thread ~{cpu} \
      --html ~{samplename}_fastplong.html --json ~{samplename}_fastplong.json \
      ~{if defined(fastplong_args) then "~{fastplong_args}" else ""}
  >>>
  output {
    File read1_trimmed = "~{samplename}_trim.fastq.gz"
    File fastplong_stats_html = "~{samplename}_fastplong.html"
    File fastplong_stats_json = "~{samplename}_fastplong.json"
    String fastplong_version = read_string("VERSION")
    String fastplong_docker = "~{docker}"
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