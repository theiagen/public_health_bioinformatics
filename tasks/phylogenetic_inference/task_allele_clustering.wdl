version 1.0

task allele_clustering {
  meta {
  description: "PulseNet 2.0 Allele Clustering algorithm"
  }
  input {
    Array[File] allele_jsons

    String tree_building_algorithm
    String distance_algorithm

    String tree_name

    Int cpu = 2
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/allele-clustering:1.0.0"
    Int memory = 4
  }
  command <<<
    # concatenate all GZIPPED JSONs into NDJSON format
    touch ~{tree_name}_concatenated_profiles.ndjson.gz
    file_array=(~{sep=' ' allele_jsons})
    for index in ${!file_array[@]}; do
      # in the case where the input files do not end in newline characters,
      # AlleleClustering.py gets angry since it expects NDJSON
      # also you cannot add newlines when concatenating gzipped files
      # so we need to unzip before concatenating within an echo for newline addition
      gunzip ${file_array[$index]} -c > temporary.json
      echo "$(cat temporary.json)" >> ~{tree_name}_concatenated_profiles.ndjson
    done

    # run AlleleClustering script
    python3 /data/AlleleClustering.py ~{tree_name}_concatenated_profiles.ndjson \
      --algorithm ~{tree_building_algorithm} \
      --distance ~{distance_algorithm} \
      --output ~{tree_name}
  >>>
  output {
    File concatenated_jsons = "~{tree_name}_concatenated_profiles.ndjson"
    File tree = "~{tree_name}.nwk"
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 1
  }
}
