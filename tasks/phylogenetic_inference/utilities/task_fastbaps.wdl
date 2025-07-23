version 1.0

task fastbaps {
  input {
    File fasta
    File tree # must be rooted
    Int? levels
    String? character # symmetric, baps, optimise.symmetric, optimise.baps
    
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/fastbaps:1.0.8-dev"
    Int disk_size = 10
    Int memory = 8
    Int cpu = 2
  }
  command <<<
    # set -euo pipefail to avoid silent failure
    set -euo pipefail

    # grab the fastbaps version
    echo $FASTBAPS_VER | tee VERSION

    # run fastbaps
    run_fastbaps \
      --input ~{fasta} \
      --phylogeny ~{tree} \
      --threads ~{cpu} \
      ~{if defined(levels) then "--levels ~{levels}" else ""} \
      ~{if defined(character) then "--character ~{character}" else ""}
  >>>
  runtime {
    docker: docker 
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
  output {
    String fastbaps_version = read_string("VERSION")
    File? fastbaps_clusters = "fastbaps_clusters.csv"
  }
}
