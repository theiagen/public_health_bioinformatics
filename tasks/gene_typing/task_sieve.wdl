version 1.0

task sieve {
  input {
    String samplename
    String plugin
    File? assembly
    File? read1
    File? read2

    File? database # tar.gz compressed database
    String? engine # blast / kma

    String? serogroup_key # key for parsing serogrouping from json output
    String? genes_key # key for parsing genes present from json output
    String? notes_key # key for parsing notes from json

    Float? min_identity
    Float? min_coverage

    String? parameters # comma-delimited list of parameters to populate

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/sieve:v0.1.0"
    Int disk_size = 50
    Int cpu = 1
    Int memory = 4
  }
  command <<<
    set -euo pipefail

    # capture version and available plugins for provenance/logs
    sieve version | tee VERSION
    sieve list-plugins

    # select input(s): assembly OR reads (single- or paired-end)
    # --input is the assembly or R1; --input2 is the R2 mate
    # reads require the KMA engine; assemblies default to BLAST
    engine="~{engine}"
    engine_args=()
    if [[ -n "~{assembly}" ]]; then
      input_args=(--input "~{assembly}")
      [[ -n "${engine}" ]] && engine_args=(--engine "${engine}")
    elif [[ -n "~{read1}" ]]; then
      input_args=(--input "~{read1}")
      if [[ -n "~{read2}" ]]; then
        input_args+=(--input2 "~{read2}")
      fi
      # enforce KMA for read inputs; error on any conflicting engine override
      if [[ -n "${engine}" && "${engine,,}" != "kma" ]]; then
        echo "ERROR: read inputs require the KMA engine, but engine='${engine}' was requested." >&2
        exit 1
      fi
      engine_args=(--engine kma)
    else
      echo "ERROR: must provide either 'assembly' or 'read1'." >&2
      exit 1
    fi

    # optional database (tar.gz of a directory holding sequences.fasta)
    db_args=()
    if [[ -n "~{database}" ]]; then
      mkdir -p sieve_db
      tar -xzvf "~{database}" -C sieve_db
      # point --db at the directory that actually contains sequences.fasta
      seqfile=$(find sieve_db -name 'sequences.fasta' | head -n1)
      if [[ -n "${seqfile}" ]]; then
        db_args=(--db "$(dirname "${seqfile}")")
      else
        db_args=(--db sieve_db)
      fi
    fi

    # optional aligner params: comma-delimited KEY=VALUE list -> repeated --param
    param_args=()
    if [[ -n "~{parameters}" ]]; then
      IFS=',' read -ra _raw_params <<< "~{parameters}"
      for p in "${_raw_params[@]}"; do
        # trim surrounding whitespace
        p="${p#"${p%%[![:space:]]*}"}"
        p="${p%"${p##*[![:space:]]}"}"
        [[ -n "${p}" ]] && param_args+=(--param "${p}")
      done
    fi

    # run sieve
    sieve screen \
      "${input_args[@]}" \
      --plugin ~{plugin} \
      --name ~{samplename} \
      "${engine_args[@]}" \
      ~{"--min-identity " + min_identity} \
      ~{"--min-coverage " + min_coverage} \
      "${db_args[@]}" \
      "${param_args[@]}" \
      --out-dir sieve/

    # always create the parsed outputs so read_string never fails
    : > SEROGROUP
    : > GENES_PRESENT
    : > NOTES

    # pass values through the environment so quotes/whitespace cannot break the python literal
    export SIEVE_JSON="sieve/~{samplename}.json"
    export SIEVE_SAMPLENAME="~{samplename}"
    export SEROGROUP_KEY="~{default='' serogroup_key}"
    export GENES_KEY="~{default='' genes_key}"
    export NOTES_KEY="~{default='' notes_key}"

    python3 <<'CODE'
    import json
    import os
    import sys

    with open(os.environ["SIEVE_JSON"]) as handle:
        report = json.load(handle)

    # all reportable fields live in the "result" object; fall back to the top level if a
    # plugin ever emits a flat report
    result = report.get("result")
    if not isinstance(result, dict):
        print("WARNING: no 'result' object in the sieve JSON; parsing top-level keys instead.", file=sys.stderr)
        result = report

    def format_value(value):
        # lists (e.g. genes_present) are reported as a comma-delimited string
        if value is None:
            return ""
        if isinstance(value, list):
            return ",".join(str(item) for item in value)
        if isinstance(value, dict):
            return json.dumps(value)
        return str(value)

    for env_var, out_file in (
        ("SEROGROUP_KEY", "SEROGROUP"),
        ("GENES_KEY", "GENES_PRESENT"),
        ("NOTES_KEY", "NOTES"),
    ):
        key = os.environ.get(env_var, "").strip()
        if not key:
            print(f"INFO: no {env_var} provided; leaving {out_file} empty.", file=sys.stderr)
            continue
        if key not in result:
            print(f"WARNING: key '{key}' not found in the sieve result for "
                  f"{os.environ['SIEVE_SAMPLENAME']}; leaving {out_file} empty.", file=sys.stderr)
            continue
        with open(out_file, "w") as out:
            out.write(format_value(result[key]) + "\n")
    CODE
  >>>
  output {
    File sieve_results = "sieve/~{samplename}.tsv"
    String sieve_serogroup = read_string("SEROGROUP")
    String sieve_genes_present = read_string("GENES_PRESENT")
    String sieve_notes = read_string("NOTES")
    String sieve_version = read_string("VERSION")
    String sieve_plugin = plugin
    String sieve_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 1 # does not take long (usually <3 min) to run stxtyper on 1 genome, preemptible is fine
    maxRetries: 3
  }
}
