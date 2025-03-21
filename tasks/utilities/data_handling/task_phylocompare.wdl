version 1.0

task phylovalidate {
  input {
    File tree1_path
    File tree2_path

    String? outgroups
    Boolean? midpoint = false
    Float? max_distance
    
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/theiaphylo:0.1.0"  # update!!!
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

    # grab the phylocompare version
    phylocompare.py --version | tee VERSION

    # set bash variables to check them for population in conditionals
    outgroups=~{outgroups}
    max_distance=~{max_distance}

    # root if outgroups are provided
    if [[ -n ${outgroups} ]]; then
      phylocompare.py ~{tree1_path} \
        ~{tree2_path} \
        --outgroup ~{outgroups} \
        --debug
    # root at the midpoint
    elif ~{midpoint}; then
      phylocompare.py ~{tree1_path} \
        ~{tree2_path} \
        --midpoint \
        --debug
    # root is not specified, allow phylocompare.py to determine
    else
      phylocompare.py ~{tree1_path} \
        ~{tree2_path} \
        --debug
    fi

    # extract the distance
    tail -1 phylo_distances.txt | cut -f 2 | tr -d ' ' > PHYLOCOMPARE_DISTANCE

    # run the comparison
    if [[ -z ${max_distance} ]]; then
      echo "NA" > PHYLOVALIDATE
    else
      python3 <<CODE
    try:
      # check if the distance is greater than the max distance
      with open('PHYLOCOMPARE_DISTANCE', 'r') as f:
        observed_distance = float(f.read().strip())
      if observed_distance > ~{max_distance}:
        with open('PHYLOVALIDATE', 'w') as out:
          out.write('FAIL')
      else:
        with open('PHYLOVALIDATE', 'w') as out:
          out.write('PASS')
    # indicates that the distance is not a float, likely a None
    except ValueError:
      with open('PHYLOVALIDATE', 'w') as out:
        out.write('FAIL')
    CODE
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
    String phylocompare_version = read_string("VERSION")
    File summary_report = "phylo_distances.txt"
    String phylo_distance = read_float("PHYLOCOMPARE_DISTANCE")
    String validation = read_string("PHYLOVALIDATE")
  }
}
