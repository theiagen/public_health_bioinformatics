version 1.0

task augur_export {
  input {
    File tree
    File? metadata
    Array[File] node_data_jsons
    String build_name

    File? auspice_config # auspice configuration file
    String? title # title to be displayed by Auspice
    File? description_md # markdown file with description of build and/or acknowledgements
    File? colors_tsv # custom color definitions, one per line
    File? lat_longs_tsv # latitudes and longitudes for geography traits
    Boolean include_root_sequence = false # export an additional json containing the root sequence used to identify mutations
  
    Int disk_size = 100
    Int memory = 64
    Int cpu = 4
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/augur:31.5.0"
  }
  command <<<
    # fail hard
    set -euo pipefail

    # prepare node_data argument
    node_data=~{sep=' ' node_data_jsons}
    if [ -z $node_data ]; then
      node_data_arg="--node-data "${node_data}
    else
      node_data_arg=""
    fi

    # check if the bash variable is empty, then check if the file is (allows for empty file inputs)
    metadata=~{metadata}
    metadata_arg=""
    if [[ ! -z $metadata ]]; then
      if [ -s $metadata ]; then
        metadata_arg="--metadata "~{metadata}
      fi
    fi

    augur export v2 \
      --tree ~{tree} \
      $metadata_arg \
      $node_data_arg \
      --output ~{build_name}_auspice.json \
      ~{"--auspice-config " + auspice_config} \
      ~{"--title " + title} \
      ~{"--description " + description_md} \
      ~{"--colors " + colors_tsv} \
      ~{"--lat-longs " + lat_longs_tsv} \
      ~{true="--include-root-sequence " false=""  include_root_sequence}
  >>>
  output {
    File auspice_json = "~{build_name}_auspice.json"
    File? root_sequence_json = "~{build_name}_auspice_root-sequence.json"
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    dx_instance_type: "mem3_ssd1_v2_x4"
    preemptible: 0
    maxRetries: 3
  }
}