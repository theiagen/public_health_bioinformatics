version 1.0

task contaminant_check {
  input {
    String expected_sequences # comma-delimited list of expected sequences, OR a key into expected_sequences_json when that is provided
    File? expected_sequences_json # optional JSON mapping of {"<NAME>": ["<SEQ1>", "<SEQ2>", ...]}; when provided, expected_sequences is used as the key to look up the list of expected sequences
    File coverage_by_sequence_json # task_mapping_stats output: coverage_by_sequence_json
    File depth_by_sequence_json # task_mapping_stats output: depth_by_sequence_json
    File reads_by_sequence_json # task_mapping_stats output: reads_by_sequence_json
    File contaminant_fasta # FASTA of contaminant sequences
    Float min_percent_coverage = 0
    Int min_depth = 0
    Int min_reads_mapped = 0

    Int? min_expected_seq # default is defined in python as all expected_sequences
    Int max_unexpected_seq = 0 # default is 0

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/pysam:1.23.1"
    Int memory = 4
    Int cpu = 2
    Int disk_size = 100
  }
  command <<<
  set -euo pipefail

  # Resolve the list of expected contaminant sequences.
  # If expected_sequences_json is provided, treat expected_sequences as a key into
  # that JSON mapping (key -> list of sequence names) and resolve it to a
  # comma-delimited list. Otherwise, use expected_sequences verbatim.
  python3 <<CODE
  import json

  expected = """~{expected_sequences}"""
  json_path = """~{default="" expected_sequences_json}"""

  if json_path:
    with open(json_path) as infile:
      mapping = {k.strip().lower(): v for k, v in json.load(infile).items()}
    # enable comma-delimited input of expected contaminant sets
    expected_list = expected.replace(" ", "").split(",")
    resolved_list = []
    for exp in expected_list:
      exp_clean = exp.strip().lower()
      if exp_clean not in mapping:
        raise KeyError(f"Key '{exp} ({exp_clean})'' not found in expected_sequences_json mapping")
      value = mapping[exp_clean]
      # if it's a json list, append as a comma-delimited string
      if isinstance(value, list):
        temp_resolved = ",".join(str(sequence) for sequence in value)
      # otherwise coerce into a compatible comma-delimited string by removing spaces
      else:
        temp_resolved = str(value).replace(" ", "")
      resolved_list.append(temp_resolved)
    # join all resolved sequences together
    resolved = ",".join(resolved_list)
  else:
    resolved = expected

  with open("EXPECTED_SEQUENCES", "w") as outfile:
    outfile.write(resolved)
  CODE

  python3 /usr/bin/contaminant_check.py \
    --expected_sequences "$(cat EXPECTED_SEQUENCES)" \
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
