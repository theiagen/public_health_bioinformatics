version 1.0

task limsport {
  input {
    File input_table
    Array[String] sample_names
    File config
    String output_delimiter = ","
    String output_prefix = "limsport"

    Boolean allow_file_parsing = false
    Int? max_file_parsing_failures

    Boolean include_qc_fails = false

    Int cpu = 1
    Int memory = 2
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/limsport:0.1.0-dev"
  }
  # limsport picks the output extension from the delimiter
  String extension = if output_delimiter == "," then "csv" else if output_delimiter == "\t" then "tsv" else "txt"
  command <<<
    set -euo pipefail

    limsport --version | cut -d' ' -f2 | tee VERSION

    limsport \
      --input ~{input_table} \
      --config ~{config} \
      --samples  ~{write_lines(sample_names)} \
      --output-delimiter "~{output_delimiter}" \
      --output-prefix "~{output_prefix}" \
      ~{true="--allow-file-parsing" false="" allow_file_parsing} \
      ~{"--max-file-parsing-failures " + max_file_parsing_failures} \
      ~{true="--include-qc-fails" false="" include_qc_fails} \
      --verbose
  >>>
  output {
    File limsport_lims_report = "~{output_prefix}_lims_report.~{extension}"
    File limsport_qc_fail_report = "~{output_prefix}_qc_fail_report.~{extension}"
    File limsport_qc_summary = "~{output_prefix}_qc_summary.~{extension}"
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
