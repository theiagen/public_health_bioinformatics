version 1.0

task mutation_context {
  input {
    File refined_tree
    File ancestral_nt_muts_json
    String build_name
    
    Int cpu = 1
    Int memory = 4
    Int disk_size = 50
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/nextstrain-mpox-mutation-context:2024-06-27"
  }
  command <<<
    # capture version information
    
    echo "DEBUG: Running mutation_context.py script now..."
    # run mutation_context.py script
    python3 /scripts/mutation_context.py \
      --tree ~{refined_tree} \
      --mutations ~{ancestral_nt_muts_json} \
      --output "~{build_name}_mpox_mutation_context.json"

    echo "DEBUG: Finished running mutation_context.py script."
  >>>
  output {
    File mutation_context_json = "~{build_name}_mpox_mutation_context.json"
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    preemptible: 1
    maxRetries: 3
  }
}