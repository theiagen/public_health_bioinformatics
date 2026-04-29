version 1.0

task contaminant_check {
  input {
    String expected_sequences # comma-delimited list of expected sequences
    File coverage_by_sequence_json # task_mapping_stats output: coverage_by_sequence_json
    File depth_by_sequence_json # task_mapping_stats output: depth_by_sequence_json
    File contaminant_fasta # FASTA of contaminant sequences
    Float min_percent_coverage = 0
    Int min_depth = 0

    Int? min_expected # default is defined in python and all expected_sequences
    Int max_unexpected = 1

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
    with open(filename, "w") as f:
      json_dump(data, f, indent=4)

  def compile_failures(passing_sequences_variable, expected_recovered_sequences, var_name):
    # sequences that passed coverage threshold and were expected
    variable_hits = passing_sequences_variable.intersection(expected_recovered_sequences)
    if variable_hits:
      print(f"DEBUG: passing {var_name}: {sorted(variable_hits)}")
    else:
      print(f"WARNING: no sequences passed {var_name} threshold")
    # sequences that failed variable threshold and were expected
    variable_missing = expected_recovered_sequences.difference(variable_hits)
    if variable_missing:
      print(f"WARNING: failing {var_name}: {sorted(variable_missing)}")
    # sequences that passed variable threshold but were not expected
    variable_extra = passing_sequences_variable.difference(expected_recovered_sequences)
    unexpected_sequences_variable = {k: variable_by_sequence[k] for k in sorted(variable_extra)}
    write_json(f"UNEXPECTED_SEQ2{var_name.upper()}", unexpected_sequences_variable)

    return variable_missing, variable_extra 
  
  # convert comma-separated string of expected sequences into a set
  expected_sequences = set([seq.strip() for seq in "~{expected_sequences}".split(",")])
  # set default to all expected_sequences
  if ~{if defined(min_expected) then "'true'" else "'false'"} == "true":
    min_expected = len(expected_sequences)
  else:
    min_expected = int(~{min_expected})

  # read in coverage and depth by sequence
  with open("~{coverage_by_sequence_json}") as f:
    coverage_by_sequence = json.load(f)
  with open("~{depth_by_sequence_json}") as f:
    depth_by_sequence = json.load(f)
  # obtain all sequences in the reference
  reference_sequences = set()
  with open("~{contaminant_fasta}") as f:
    for line in f:
      if line.startswith(">"):
        # acquire the sequence from the header by removing descriptions and the ">"
        seq = line.split()[0][1:]
        reference_sequences.add(seq)

  # check if any expected sequences are present above the specified thresholds
  passing_sequences_coverage_prep = set([seq for seq, cov in coverage_by_sequence.items() if cov >= ~{min_percent_coverage}])
  passing_sequences_depth_prep = set([seq for seq, depth in depth_by_sequence.items() if depth >= ~{min_depth}])

  failing_sequences_coverage = set([seq for seq, cov in coverage_by_sequence.items() if not cov > 0])
  failing_sequences_depth = set([seq for seq, depth in depth_by_sequence.items() if not depth > 0])
  failing_sequences = failing_sequences_coverage.union(failing_sequences_depth)

  passing_sequences_coverage = passing_sequences_coverage_prep.difference(failing_sequences_coverage)
  passing_sequences_depth = passing_sequences_depth_prep.difference(failing_sequences_depth)
  passing_sequences = passing_sequences_coverage.union(passing_sequences_depth)

  # assumes coverage_by_sequence and depth_by_sequence have same keys
  # sequences that were desired to be identified, but not recovered
  expected_unrecovered_sequences = expected_sequences.difference(reference_sequences)
  # sequences that were expected and recovered
  expected_recovered_sequences = expected_sequences.difference(expected_unrecovered_sequences)
  expected_sequences_coverage = {k: coverage_by_sequence[k] for k in sorted(expected_recovered_sequences)}
  write_json("EXPECTED_SEQ2COVERAGE.json", expected_sequences_coverage)
  expected_sequences_depth = {k: depth_by_sequence[k] for k in sorted(expected_recovered_sequences)}
  write_json("EXPECTED_SEQ2DEPTH.json", expected_sequences_depth)

  # check results
  print(f"DEBUG: expected sequences: {sorted(expected_sequences)}")
  coverage_missing, coverage_extra = compile_failures(passing_sequences_coverage, expected_recovered_sequences, "coverage")
  depth_missing, depth_extra = compile_failures(passing_sequences_depth, expected_recovered_sequences, "depth")

  # total unexpected sequences that passed thresholds
  unexpected_sequences = sorted(coverage_extra.union(depth_extra))
  if unexpected_sequences:
    print(f"WARNING: extraneous sequences detected: {sorted(unexpected_sequences)}")

  seq2fail = defaultdict(list)
  for seq in coverage_missing:
    if seq in failing_sequences or seq not in coverage_by_sequence:
      seq2fail[seq].append("not detected")
    else:
      seq2fail[seq].append(f"insufficient coverage ({coverage_by_sequence[seq]})")
  for seq in depth_missing:
    # we already accounted for the sequence not being detected outright, so this is a depth-check
    if seq not in failing_sequences and seq in depth_by_sequence:
      seq2fail[seq].append(f"insufficient depth ({depth_by_sequence[seq]})")
  # these sequences are missing from the reference because they were expected
  # but not detected/accounted for in the reference FASTA
  for seq in expected_unrecovered_sequences:
    seq2fail[seq].append("missing from reference")

  # populate a status string
  with open("STATUS", "w") as f:
    # check if a pass/fail threshold was infringed 
    if len(seq2fail) < min_expected or len(unexpected_sequences) > int(~{max_unexpected}):
      status_string = "FAIL: "
      # too few expected sequences recovered
      if len(seq2fail) < min_expected:
        for seq, fail_reasons in sorted(seq2fail.items(), key=lambda x: x[0]):
          status_string += f"{seq} - {', '.join(fail_reasons)}; "
      # too many unexpected sequences recovered
      else:
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
    Map[String, Float] unexpected_coverage_by_sequence = read_json("UNEXPECTED_SEQ2COVERAGE.json")
    Map[String, Float] unexpected_depth_by_sequence = read_json("UNEXPECTED_SEQ2DEPTH.json")
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