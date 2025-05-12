version 1.0

task root_phylo {
  input {
    File tree
    String? outgroups
    Boolean? midpoint = false
    
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/theiaphylo:0.1.1"
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
    TheiaPhylo.py --version | tee VERSION

    # set bash variable to check for population in conditionals
    outgroups=~{outgroups}

    # prepare an output file path
    rooted_tree=$(sed -Ee 's/(.*)\.[^\.]+$/\1.rooted.nwk/' <<< ${tree})

    # root if outgroups are provided
    if [[ -n ${outgroups} ]]; then
      TheiaPhylo.py ${tree} \
        --outgroup ~{outgroups} \
        --output ${rooted_tree}
    # root at the midpoint
    elif ~{midpoint}; then
      TheiaPhylo.py ${tree} \
        --midpoint \
        --output ${rooted_tree}
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
    File rooted_tree = "${rooted_tree}"
  }
}
