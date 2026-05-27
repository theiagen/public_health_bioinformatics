version 1.0

task contaminant_check {
  input {
    String expected_sequences # comma-delimited list of expected sequences
    File coverage_by_sequence_json # task_mapping_stats output: coverage_by_sequence_json
    File depth_by_sequence_json # task_mapping_stats output: depth_by_sequence_json
    File reads_by_sequence_json # task_mapping_stats output: reads_by_sequence_json
    File contaminant_fasta # FASTA of contaminant sequences
    Float min_percent_coverage = 0
    Int min_depth = 0
    Int min_reads_mapped = 0

    Int? min_expected_seq # default is defined in python as all expected_sequences
    Int max_unexpected_seq = 0 # default is 0

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/pysam:1.23"
    Int memory = 4
    Int cpu = 2
    Int disk_size = 100
  }
  command <<<
  python3 /usr/bin/contaminant_check.py \
    --expected_sequences ~{expected_sequences} \
    --coverage_by_sequence_json ~{coverage_by_sequence_json} \
    --depth_by_sequence_json ~{depth_by_sequence_json} \
    --reads_by_sequence_json ~{reads_by_sequence_json} \
    --contaminant_fasta ~{contaminant_fasta} \
    --minimum_percent_coverage ~{min_percent_coverage} \
    --minimum_depth ~{min_depth} \
    --minimum_reads_mapped ~{min_reads_mapped} \
    --maximum_unexpected_sequences ~{max_unexpected_seq} \
    ~{if defined(min_expected_seq) then "--minimum_expected_sequences ~{min_expected_seq}" else ""}
  >>>
  output {
    String contaminant_check_status = read_string("STATUS")
    Map[String, Float] expected_coverage_by_sequence = read_json("EXPECTED_SEQ2COVERAGE.json")
    Map[String, Float] expected_depth_by_sequence = read_json("EXPECTED_SEQ2DEPTH.json")
    Map[String, Float] expected_reads_by_sequence = read_json("EXPECTED_SEQ2READS.json")
    Map[String, Float] unexpected_coverage_by_sequence = read_json("UNEXPECTED_SEQ2COVERAGE.json")
    Map[String, Float] unexpected_depth_by_sequence = read_json("UNEXPECTED_SEQ2DEPTH.json")
    Map[String, Float] unexpected_reads_by_sequence = read_json("UNEXPECTED_SEQ2READS.json")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}
