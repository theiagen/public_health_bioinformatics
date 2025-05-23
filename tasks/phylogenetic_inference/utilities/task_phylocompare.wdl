version 1.0

task phylovalidate {
  input {
    File tree1
    File tree2
    Float? max_distance
    
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/theiaphylo:0.1.7-dev"
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
    phylocompare --version | tee VERSION

    # set clean tree PATHs
    tree1_clean=~{tree1}.clean
    tree2_clean=~{tree2}.clean

    # clean the trees, report if they are bifurcating
    Rscript /theiaphylo/theiaphylo/clean_phylo.R ~{tree1} > ${tree1_clean} 2> >(cut -f 2 -d ' ' > TREE1_BIFURCATING)
    Rscript /theiaphylo/theiaphylo/clean_phylo.R ~{tree2} > ${tree2_clean} 2> >(cut -f 2 -d ' ' > TREE2_BIFURCATINGH

    # set bash variables to check them for population in conditionals
    max_distance=~{max_distance}

    # run comparison
    phylocompare ${tree1_clean} ${tree2_clean} \
        --debug

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
    with open('TREE1_BIFURCATING', 'r') as f:
      tree1_bifurcating = f.read().strip()
    with open('TREE2_BIFURCATING', 'r') as f:
      tree2_bifurcating = f.read().strip()
    if tree1_bifurcating == 'FALSE' or tree2_bifurcating == 'FALSE':
      with open('PHYLOCOMPARE_FLAG', 'a') as out:
        out.write('polytomy')
    else:
      with open('PHYLOCOMPARE_FLAG', 'a') as out:
        out.write('')
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
    File tree1_clean = "~{tree1}.clean"
    File tree2_clean = "~{tree2}.clean"
    String phylo_distance = read_float("PHYLOCOMPARE_DISTANCE")
    String phylo_validation = read_string("PHYLOVALIDATE")
    String phylo_flag = read_string("PHYLOCOMPARE_FLAG")
  }
}
