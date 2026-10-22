version 1.0

task limsport {
  input {
    File input_table
    File config
    Array[String] sample_names
    String output_delimiter = ","
    Boolean allow_file_parsing = false
    Int? max_file_parsing_failures
    String output_prefix = "limsport"

    Int cpu = 1
    Int memory = 2
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/limsport:0.1.0-dev"
  }
  command <<<
    set -euo pipefail

    limsport --version | cut -d' ' -f2 | tee VERSION

    limsport \
      --input ~{input_table} \
      --config ~{config} \
      --samples  ~{write_lines(sample_names)} \
      --delimiter "~{output_delimiter}" \
      ~{true="--allow-file-parsing" false="" allow_file_parsing} \
      ~{"--max-file-parsing-failures " + max_file_parsing_failures} \
      --output "~{output_prefix}.txt" \
      --qc-report "~{output_prefix}_qc_report.tsv"
  >>>
  output {
    File limsport_table = "~{output_prefix}.txt"
    File limsport_qc_report = "~{output_prefix}_qc_report.tsv"
    String limsport_version = read_string("VERSION")
    String limsport_docker = docker
  }
  runtime {
    docker: docker
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 1
  }
}
