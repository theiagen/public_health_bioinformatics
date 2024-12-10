version 1.0

task medaka_consensus {
  input {
    File unpolished_fasta
    String samplename
    File read1
    Boolean auto_model = true       # Automatically resolve model by inspecting input
    String? medaka_model_override  # Optional user-specified model
    Int cpu = 4
    Int memory = 16
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/medaka:2.0.1"
  }
  command <<<
    set -euo pipefail

    medaka --version | tee MEDAKA_VERSION

    # Initialize selected_model variable
    selected_model=""

    # Check if auto_model is enabled
    if [[ "~{auto_model}" == "true" ]]; then
      echo "Attempting automatic model selection..."
      medaka tools resolve_model --auto_model consensus ~{read1} > auto_model.txt || true
      auto_model=$(cat auto_model.txt || echo "")

      # Alert if automatic model selection fails
      if [[ -z "$auto_model" ]]; then
        echo "Warning: Automatic model selection failed. Defaulting to fallback model."
      fi
    fi

    # Use auto_model if detected, otherwise fallback to user-provided or default
    if [[ -n "$auto_model" ]]; then
      selected_model="$auto_model"
    elif [[ -n "~{medaka_model_override}" ]]; then
      selected_model="~{medaka_model_override}"
    else
      selected_model="r1041_e82_400bps_sup_v5.0.0"
    fi

    echo "Using Medaka model for polishing: $selected_model"

    # Perform Medaka polishing
    medaka_consensus \
      -i ~{read1} \
      -d ~{unpolished_fasta} \
      -o . \
      -m $selected_model \
      -t ~{cpu}

    mv consensus.fasta ~{samplename}.polished.fasta
  >>>
  output {
    File medaka_fasta = "~{samplename}.polished.fasta"
    String medaka_version = read_string("MEDAKA_VERSION")
  }
  runtime {
    docker: "~{docker}"
    cpu: cpu
    memory: "~{memory} GB"
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    maxRetries: 1
    preemptible: 0
  }
}
