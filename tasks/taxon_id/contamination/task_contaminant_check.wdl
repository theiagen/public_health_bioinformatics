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

    Int? min_expected_seq # default is defined in python and all expected_sequences
    Int max_unexpected_seq = 1

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/pysam:1.23"
    Int memory = 4
    Int cpu = 2
    Int disk_size = 100
  }
  command <<<
  python3 <<CODE
  import json
  from collections import defaultdict

  def write_json(filename, data):
    """Write a JSON file compatible with WDL"""
    with open(filename, "w") as f:
      if data:
        json.dump(data, f, indent=4)
      else:
        # spoof Cromwell (Terra WDL)
        f.write('{"": 0}')

  def apply_thresholds(variable_by_sequence, min_value):
    """Apply thresholds for a variable while excluding empty values as missing"""
    passing_sequences_prep = set([seq for seq, val in variable_by_sequence.items() if val >= min_value])
    failing_sequences = set([seq for seq, val in variable_by_sequence.items() if not val > 0])
    passing_sequences = passing_sequences_prep.difference(failing_sequences)
    return failing_sequences, passing_sequences

  def compile_failures(passing_sequences, expected_recovered_sequences, variable_by_sequence, var_name):
    """Compile failing sequences based on those that are expected v. unexpected"""
    # sequences that passed coverage threshold and were expected
    hit_sequences = passing_sequences.intersection(expected_recovered_sequences)
    if hit_sequences:
      print(f"DEBUG: passing {var_name}: {sorted(hit_sequences)}")
    else:
      print(f"WARNING: no sequences passed {var_name} threshold")
    # sequences that failed variable threshold and were expected
    missing_sequences = expected_recovered_sequences.difference(hit_sequences)
    if missing_sequences:
      print(f"WARNING: failing {var_name}: {sorted(missing_sequences)}")
    # sequences that passed variable threshold but were not expected
    extra_sequences = passing_sequences.difference(expected_recovered_sequences)
    unexpected_sequences = {seq: variable_by_sequence[seq] for seq in sorted(extra_sequences)}
    write_json(f"UNEXPECTED_SEQ2{var_name.upper()}.json", unexpected_sequences)
    return missing_sequences, extra_sequences 

  def annotate_failures(seq2fail, variable_missing, failing_sequences, variable_by_sequence, variable):
    """Annotate failure type for the status output"""
    for seq in variable_missing:
      if seq in failing_sequences or seq not in variable_by_sequence:
        seq2fail[seq].append("not detected")
      else:
        seq2fail[seq].append(f"insufficient {variable} ({variable_by_sequence[seq]})")
    return seq2fail

  # convert comma-separated string of expected sequences into a set
  expected_sequences = set([seq.strip() for seq in "~{expected_sequences}".split(",") if seq.strip()])
  # set default to all expected_sequences
  if ~{if defined(min_expected_seq) then "True" else "False"}:
    min_expected_seq = int(~{min_expected_seq})
    if min_expected_seq > len(expected_sequences):
      print(f"ERROR: min_expected_seq ({min_expected_seq}) exceeds number of expected_sequences ({len(expected_sequences)}); setting min_expected_seq to {len(expected_sequences)}")
      min_expected_seq = len(expected_sequences)
  else:
    min_expected_seq = len(expected_sequences)
  print(f"DEBUG: expecting minimum {min_expected_seq} sequences")

  # read in coverage and depth by sequence
  with open("~{coverage_by_sequence_json}") as f:
    coverage_by_sequence = json.load(f)
  with open("~{depth_by_sequence_json}") as f:
    depth_by_sequence = json.load(f)
  with open("~{reads_by_sequence_json}") as f:
    reads_by_sequence = json.load(f)
  # obtain all sequences in the reference
  reference_sequences = set()
  with open("~{contaminant_fasta}") as f:
    for line in f:
      if line.startswith(">"):
        # acquire the sequence from the header by removing descriptions and the ">"
        seq = line.split()[0][1:]
        reference_sequences.add(seq)

  # check if any expected sequences are present above the specified thresholds
  failing_sequences_coverage, passing_sequences_coverage = apply_thresholds(coverage_by_sequence, ~{min_percent_coverage})
  failing_sequences_depth, passing_sequences_depth = apply_thresholds(depth_by_sequence, ~{min_depth})
  failing_sequences_reads, passing_sequences_reads = apply_thresholds(reads_by_sequence, ~{min_reads_mapped})

  failing_sequences = failing_sequences_coverage.union(failing_sequences_depth).union(failing_sequences_reads)
  passing_sequences = passing_sequences_coverage.union(passing_sequences_depth).union(passing_sequences_reads)

  # sequences that were desired to be identified, but not recovered
  expected_unrecovered_sequences = expected_sequences.difference(reference_sequences)
  # sequences that were expected and recovered
  expected_recovered_sequences = expected_sequences.difference(expected_unrecovered_sequences)

  # write outputs for recovered expected sequences 
  expected_sequences_coverage = {}
  expected_sequences_depth = {}
  expected_sequences_reads = {}
  for seq in sorted(expected_recovered_sequences):
    expected_sequences_coverage[seq] = coverage_by_sequence[seq]
    expected_sequences_depth[seq] = depth_by_sequence[seq]
    expected_sequences_reads[seq] = reads_by_sequence[seq]
  write_json("EXPECTED_SEQ2COVERAGE.json", expected_sequences_coverage)
  write_json("EXPECTED_SEQ2DEPTH.json", expected_sequences_depth)
  write_json("EXPECTED_SEQ2READS.json", expected_sequences_reads)

  # check results and write outputs for unexpected sequences
  print(f"DEBUG: expected sequences: {sorted(expected_sequences)}")
  coverage_missing, coverage_extra = compile_failures(passing_sequences_coverage, 
                                                      expected_recovered_sequences, 
                                                      coverage_by_sequence, 
                                                      "coverage")
  depth_missing, depth_extra = compile_failures(passing_sequences_depth, 
                                                expected_recovered_sequences, 
                                                depth_by_sequence, 
                                                "depth")
  reads_missing, reads_extra = compile_failures(passing_sequences_reads, 
                                                expected_recovered_sequences, 
                                                reads_by_sequence, 
                                                "reads")

  # total unexpected sequences that passed thresholds
  unexpected_sequences = sorted(coverage_extra.union(depth_extra).union(reads_extra))
  if unexpected_sequences:
    print(f"WARNING: extraneous sequences detected: {sorted(unexpected_sequences)}")

  # annotate failing sequences
  seq2fail = defaultdict(list)
  seq2fail = annotate_failures(seq2fail, coverage_missing, failing_sequences, coverage_by_sequence, "coverage")
  seq2fail = annotate_failures(seq2fail, depth_missing, failing_sequences, depth_by_sequence, "depth")
  seq2fail = annotate_failures(seq2fail, reads_missing, failing_sequences, reads_by_sequence, "reads mapped")
  # these sequences are missing from the reference because they were expected
  # but not detected/accounted for in the reference FASTA
  for seq in expected_unrecovered_sequences:
    seq2fail[seq].append("missing from reference")

  # populate a status string
  with open("STATUS", "w") as f:
    # check if a pass/fail threshold was infringed 
    if len(expected_recovered_sequences) < min_expected_seq or len(unexpected_sequences) > int(~{max_unexpected_seq}):
      status_string = "FAIL: "
      # too few expected sequences recovered
      if len(expected_recovered_sequences) < min_expected_seq:
        for seq, fail_reasons in sorted(seq2fail.items(), key=lambda x: x[0]):
          status_string += f"{seq} - {', '.join(fail_reasons)}; "
      # too many unexpected sequences recovered
      if len(unexpected_sequences) > int(~{max_unexpected_seq}):
        for seq in unexpected_sequences:
          status_string += f"{seq} - extra sequence; "
      status_string = status_string.strip("; ")
      f.write(status_string)
    else:
      f.write("PASS")
  CODE
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