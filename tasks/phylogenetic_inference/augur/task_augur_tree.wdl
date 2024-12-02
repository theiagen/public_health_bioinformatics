version 1.0

task augur_tree {
  input {
    File aligned_fasta
    String build_name
    String method = "iqtree" # possible choices: fasttree, raxml, iqtree
    String substitution_model = "GTR" # only available for iqtree
    File? exclude_sites # file name of one-based sites to exclude for raw tree building
    String? tree_builder_args # additional tree builder arguments
    Boolean override_default_args = false # override default tree builder args instead of augmenting them

    Int cpu = 64
    Int memory = 32
    Int disk_size = 750
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/augur:22.0.2--pyhdfd78af_0"
  }
  command <<<
    set -euo pipefail
    
    # capture version information
    augur version > VERSION
    echo
    echo "mafft version:"
    mafft --version 2>&1 | tee MAFFT_VERSION
    echo
    echo "iqtree version:"
    iqtree --version | grep version | sed 's/.*version/version/;s/ for Linux.*//' | tee IQTREE_VERSION
    echo
    echo "Running augur tree now..."

    AUGUR_RECURSION_LIMIT=10000 augur tree \
      --alignment "~{aligned_fasta}" \
      --output "~{build_name}_~{method}.nwk" \
      --method "~{method}" \
      --substitution-model ~{substitution_model} \
      ~{"--exclude-sites " + exclude_sites} \
      ~{"--tree-builder-args " + tree_builder_args} \
      ~{true="--override-default-args" false="" override_default_args} \
      --nthreads auto

    # If iqtree, get the model used
    if [ "~{method}" == "iqtree" ]; then
      if [ "~{substitution_model}" == "auto" ]; then
        FASTA_BASENAME=$(basename ~{aligned_fasta} .fasta)
        FASTA_DIR=$(dirname ~{aligned_fasta})
        MODEL=$(grep "Best-fit model:" ${FASTA_DIR}/*${FASTA_BASENAME}-delim.iqtree.log | sed 's|Best-fit model: ||g;s|chosen.*||' | tr -d '\n\r')
      else
        MODEL="~{substitution_model}"
      fi
      echo "$MODEL" > FINAL_MODEL.txt
    else
      echo "" > FINAL_MODEL.txt
    fi

    echo 
    echo "DEBUG: FINAL_MODEL.txt is: $(cat FINAL_MODEL.txt)"
  >>>

  output {
    File aligned_tree  = "~{build_name}_~{method}.nwk"
    String augur_version = read_string("VERSION")
    String mafft_version = read_string("MAFFT_VERSION")
    String iqtree_version = read_string("IQTREE_VERSION")
    String iqtree_model_used = read_string("FINAL_MODEL.txt")
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x36"
    preemptible: 0
    maxRetries: 3
  }
}