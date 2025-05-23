version 1.0

task root_phylo {
  input {
    File tree
    String? outgroups
    Boolean? midpoint = false
    
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/theiaphylo:0.1.6"
    Int disk_size = 10
    Int memory = 4
    Int cpu = 1
  }
  meta {
    # added so that call caching is always turned off
    volatile: true
  }
  command <<<
    # set -euo pipefail to avoid silent failure
    set -euo pipefail

    # grab the TheiaPhylo version
    theiaphylo --version | tee VERSION

    # set bash variable to check for population in conditionals
    outgroups=~{outgroups}

    # prepare an output file path
    rooted_tree=~{tree}.rooted.nwk

    # root if outgroups are provided
    if [[ -n ${outgroups} ]]; then
      theiaphylo ${tree} \
        --outgroup ~{outgroups} \
        --output ${rooted_tree}
    # root at the midpoint
    elif ~{midpoint}; then
      theiaphylo ${tree} \
        --midpoint \
        --output ${rooted_tree}
    fi
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
    String theiaphylo_version = read_string("VERSION")
    File rooted_tree = "~{tree}.rooted.nwk"
  }
}
