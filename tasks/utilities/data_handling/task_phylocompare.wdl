version 1.0

task phylovalidate {
  input {
    File tree1_path
    File tree2_path

    String? root_tips
    Boolean? unrooted = false
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

    # grab the ete3 version
    phylocompare.py --version | tee VERSION

    # set bash variables to check them for population in conditionals
    root_tips=~{root_tips}
    max_distance=~{max_distance}

    # root if outgroups are provided
    if [[ -n ${root_tips} ]]; then
      phylocompare.py ~{tree1_path} ~{tree2_path} \
        --outgroup ~{root_tips} \
        --debug
    # root at the midpoint
    elif ~{midpoint}; then
      phylocompare.py ~{tree1_path} ~{tree2_path} \
        --midpoint \
        --debug
    # append unrooted if provided, otherwise assume tree is prerooted
    else
      phylocompare.py ~{tree1_path} ~{tree2_path} \
        ~{true="--unrooted" false="" unrooted} \
        --debug
    fi

    # extract the distance
    tail -1 phylo_distances.txt | cut -f 2 | tr -d ' ' > PHYLOCOMPARE_DISTANCE

    # run the comparison
    if [[ -z ${max_distance} ]]; then
      echo "NA" > phylovalidate
    else
      python3 -c "if float(open('PHYLOCOMPARE_DISTANCE', 'r').read().strip()) > ~{max_distance}: open('phylovalidate', 'w').write('FAIL'); else: open('phylovalidate', 'w').write('PASS')"
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
    Float phylo_distance = read_float("PHYLOCOMPARE_DISTANCE")
    String validation = read_string("phylovalidate")
  }
}
