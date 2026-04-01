version 1.0

task contaminant_check {
  input {
    String expected_sequences # comma-delimited list of expected sequences
    File coverage_by_sequence_json # task_mapping_stats output: coverage_by_sequence_json
    File depth_by_sequence_json # task_mapping_stats output: depth_by_sequence_json
    File? cov_stats # task_mapping_stats output: cov_stats_txt
    Float min_percent_coverage = 0
    Int min_depth = 0

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/pysam:1.23"
    Int memory = 4
    Int cpu = 2
    Int disk_size = 100
  }
  command <<<
  python3 <<CODE
  import json
  from collections import defaultdict

  # convert comma-separated string of expected sequences into a set
  expected_sequences = set([seq.strip() for seq in "~{expected_sequences}".split(",")])

  # read in coverage and depth by sequence
  with open("~{coverage_by_sequence_json}") as f:
    coverage_by_sequence = json.load(f)
  with open("~{depth_by_sequence_json}") as f:
    depth_by_sequence = json.load(f)
  if "~{if defined(cov_stats) then 'true' else 'false'}" == "true":
    with open("~{cov_stats}") as f:
      all_sequences = set([line.split("\t")[0].strip() for line in f.readlines()[1:]])
      if "" in all_sequences:
        all_sequences.remove("")
  else:
    all_sequences = set()

  # check if any expected sequences are present above the specified thresholds
  detected_sequences_coverage_prep = set([seq for seq, cov in coverage_by_sequence.items() if cov >= ~{min_percent_coverage}])
  detected_sequences_depth_prep = set([seq for seq, depth in depth_by_sequence.items() if depth >= ~{min_depth}])

  undetected_sequences_coverage = set([seq for seq, cov in coverage_by_sequence.items() if not cov > 0])
  undetected_sequences_depth = set([seq for seq, depth in depth_by_sequence.items() if not depth > 0])
  undetected_sequences = undetected_sequences_coverage.union(undetected_sequences_depth)
  # assumes coverage_by_sequence and depth_by_sequence have same keys
  missing_from_input = expected_sequences.difference(all_sequences)

  detected_sequences_coverage = detected_sequences_coverage_prep.difference(undetected_sequences_coverage)
  detected_sequences_depth = detected_sequences_depth_prep.difference(undetected_sequences_depth)

  expected_id_sequences = expected_sequences.difference(missing_from_input)

  # check results
  print(f"DEBUG: expected sequences: {sorted(expected_sequences)}")
  coverage_hits = detected_sequences_coverage.intersection(expected_id_sequences)
  if coverage_hits:
    print(f"DEBUG: passing coverage: {sorted(coverage_hits)}")
  else:
    print("WARNING: no sequences passed coverage threshold")
  coverage_missing = expected_id_sequences.difference(coverage_hits)
  if coverage_missing:
    print(f"WARNING: failing coverage: {sorted(coverage_missing)}")
  coverage_extra = detected_sequences_coverage.difference(expected_id_sequences)
  depth_hits = detected_sequences_depth.intersection(expected_id_sequences)
  if depth_hits:
    print(f"DEBUG: passing depth: {sorted(depth_hits)}")
  else:
    print("WARNING: no sequences passed depth threshold")
  depth_missing = expected_id_sequences.difference(depth_hits)
  if depth_missing:
    print(f"WARNING: failing depth: {sorted(depth_missing)}")
  depth_extra = detected_sequences_depth.difference(expected_id_sequences)
  extraneous_sequences = coverage_extra.union(depth_extra)
  if extraneous_sequences:
    print(f"WARNING: extraneous sequences detected: {sorted(extraneous_sequences)}")

  seq2fail = defaultdict(list)
  for seq in coverage_missing:
    if seq in undetected_sequences or seq not in coverage_by_sequence:
      seq2fail[seq].append("not detected")
    else:
      seq2fail[seq].append(f"insufficient coverage ({coverage_by_sequence[seq]})")
  for seq in depth_missing:
    if seq not in undetected_sequences:
      seq2fail[seq].append(f"insufficient depth ({depth_by_sequence[seq]})")
  for seq in extraneous_sequences:
    seq2fail[seq].append("extra sequence")
  # identify missing from reference if the capacity exists
  if all_sequences:
    for seq in missing_from_input:
      seq2fail[seq].append("missing from reference")

  # populate a status string
  with open("STATUS", "w") as f:
    if seq2fail:
      status_string = "FAIL: "
      for seq, fail_reasons in sorted(seq2fail.items(), key=lambda x: x[0]):
        status_string += f"{seq} - {', '.join(fail_reasons)}; "
      status_string = status_string.strip("; ")
      f.write(status_string)
    else:
      f.write("PASS")
  CODE
  >>>
  output {
    String contaminant_check_status = read_string("STATUS")
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