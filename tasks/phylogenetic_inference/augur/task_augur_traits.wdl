version 1.0

task augur_traits {
  input {
    File refined_tree
    File? metadata
    File? weights
    #Boolean confidence = true
    String? metadata_id_columns
    String columns
    String build_name

    Int memory = 30
    Int cpu = 4
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/augur:31.5.0"
  }
  command <<<
    AUGUR_RECURSION_LIMIT=10000 augur traits \
      --tree "~{refined_tree}" \
      ~{'--metadata ' + metadata} \
      ~{'--columns {' + columns + '}'} \
      --confidence \
      ~{'--metadata-id-columns ' + metadata_id_columns} \
      ~{'--weights ' + weights} \
      --output-node-data "~{build_name}_traits.json"
  >>>
  output {
    File traits_assignments_json = "~{build_name}_traits.json"
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" 
    dx_instance_type: "mem3_ssd2_x4"
    preemptible: 0
    maxRetries: 3
  }
}