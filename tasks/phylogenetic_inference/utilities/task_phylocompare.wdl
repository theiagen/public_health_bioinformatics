version 1.0

task phylovalidate {
  input {
    File tree1
    File tree2
    Float? max_distance
    
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/theiaphylo:0.1.7"
    Int disk_size = 10
    Int memory = 4
    Int cpu = 1
  }
  String tree1_cleaned = basename(tree1) + ".clean"
  String tree2_cleaned = basename(tree2) + ".clean"
  command <<<
    # set -euo pipefail to avoid silent failure
    set -euo pipefail

    # grab the phylocompare version
    phylocompare --version | tee VERSION

    # clean the trees, report if they are bifurcating
    Rscript /theiaphylo/theiaphylo/clean_phylo.R ~{tree1} > ~{tree1_cleaned} 2> >(cut -f 2 -d ' ' > TREE1_BIFURCATING)
    Rscript /theiaphylo/theiaphylo/clean_phylo.R ~{tree2} > ~{tree2_cleaned} 2> >(cut -f 2 -d ' ' > TREE2_BIFURCATING)

    # set bash variables to check them for population in conditionals
    max_distance=~{max_distance}

    # run comparison
    phylocompare ~{tree1_cleaned} ~{tree2_cleaned} \
        --debug \
        2> PHYLOCOMPARE_STDERR

    # extract errors
    grep -Po "ERROR.*" PHYLOCOMPARE_STDERR > PHYLOCOMPARE_ERRORS

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
      with open('PHYLOCOMPAREDISTANCE', 'w') as out:
        out.write('None')
    with open('TREE1_BIFURCATING', 'r') as f:
      tree1_bifurcating = f.read().strip()
    with open('TREE2_BIFURCATING', 'r') as f:
      tree2_bifurcating = f.read().strip()

    phylocompare_flags = []
    if tree1_bifurcating == 'FALSE' or tree2_bifurcating == 'FALSE':
      phylocompare_flags.append('polytomy')
    with open('PHYLOCOMPARE_ERRORS', 'r') as f:
      errors = set(x.strip() for x in f)
      if "ERROR - Error comparing trees: number of edges must be equal" in errors:
        phylocompare_flags.append('edge_count_mismatch')
    with open('PHYLOCOMPARE_FLAG', 'w') as out:
      out.write(', '.join(phylocompare_flags))
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
    File tree1_clean = "~{tree1_cleaned}"
    File tree2_clean = "~{tree2_cleaned}"
    String phylovalidate_distance = read_string("PHYLOCOMPARE_DISTANCE")
    String phylovalidate_validation = read_string("PHYLOVALIDATE")
    String phylovalidate_flag = read_string("PHYLOCOMPARE_FLAG")
  }
}
