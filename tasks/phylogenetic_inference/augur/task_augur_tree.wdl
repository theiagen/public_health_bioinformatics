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

    Int cpus = 64
    Int disk_size = 750
  }
  command <<<
    AUGUR_RECURSION_LIMIT=10000 augur tree \
      --alignment "~{aligned_fasta}" \
      --output "~{build_name}_~{method}.nwk" \
      --method "~{method}" \
      --substitution-model ~{substitution_model} \
      ~{"--exclude-sites " + exclude_sites} \
      ~{"--tree-builder-args " + tree_builder_args} \
      ~{true="--override-default-args" false="" override_default_args} \
      --nthreads auto
  >>>
  output {
    File aligned_tree  = "~{build_name}_~{method}.nwk"
  }
  runtime {
    docker: "quay.io/biocontainers/augur:22.0.2--pyhdfd78af_0"
    memory: "32 GB"
    cpu: cpus
    disks: "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x36"
    preemptible: 0
    maxRetries: 3
  }
}