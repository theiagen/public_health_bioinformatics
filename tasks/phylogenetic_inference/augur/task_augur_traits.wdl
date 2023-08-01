version 1.0

task augur_traits {
  input {
    File refined_tree
    File metadata
    File? weights
    #Boolean confidence = true
    String? metadata_id_columns
    String columns
    String build_name

    Int mem_size = 30
    Int disk_size = 100
  }
  command <<<
    AUGUR_RECURSION_LIMIT=10000 augur traits \
      --tree "~{refined_tree}" \
      --metadata "~{metadata}" \
      --columns "~{columns}" \
      --confidence \
      ~{'--metadata-id-columns ' + metadata_id_columns} \
      ~{'--weights ' + weights} \
      --output-node-data "~{build_name}_traits.json"
  >>>
  output {
    File traits_assignments_json = "~{build_name}_traits.json"
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/biocontainers/augur:22.0.2--pyhdfd78af_0"
    memory: mem_size + " GB"
    cpu: 4
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" 
    dx_instance_type: "mem3_ssd2_x4"
    preemptible: 0
    maxRetries: 3
  }
}