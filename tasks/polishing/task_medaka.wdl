version 1.0

task medaka {
  input {
    File unpolished_fasta
    String samplename
    File read1
    Boolean auto_model = true # Enable automatic Medaka model selection
    String medaka_model = "r1041_e82_400bps_sup_v5.0.0" # Default model if auto_model is disabled or no model is resolved

    Int cpu = 4
    Int memory = 16
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/medaka:2.0.1"
  }
  command <<<    
    set -euo pipefail

    medaka --version | tee MEDAKA_VERSION

    # Attempt automatic model resolution if enabled
    if [[ "~{auto_model}" == "true" ]]; then
      resolved_model=$(medaka tools resolve_model --auto_model consensus ~{read1} 2>/dev/null || echo "~{medaka_model}")
      echo "Using model: $resolved_model"
    else
      resolved_model="~{medaka_model}"
      echo "Using specified model: $resolved_model"
    fi

    echo "$resolved_model" > MEDAKA_MODEL

    # Perform Medaka polishing
    medaka_consensus \
      -i "~{read1}" \
      -d "~{unpolished_fasta}" \
      -o . \
      -m "$resolved_model" \
      -t "~{cpu}"

    mv consensus.fasta ~{samplename}.polished.fasta
  >>>
  output {
    File medaka_fasta = "~{samplename}.polished.fasta"
    String medaka_version = read_string("MEDAKA_VERSION")
    String resolved_medaka_model = read_string("MEDAKA_MODEL")
  }
  runtime {
    docker: "~{docker}"
    cpu: cpu
    memory: "~{memory} GB"
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}
