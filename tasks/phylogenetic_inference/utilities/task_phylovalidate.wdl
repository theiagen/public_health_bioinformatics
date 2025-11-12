version 1.0

task phylovalidate {
  input {
    File tree1
    File tree2
    Float? max_distance
    Boolean resolve_tip_discrepancies = true
    
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/theiaphylo:0.2.0"
    Int disk_size = 10
    Int memory = 4
    Int cpu = 1
  }
  String tree1_cleaned = basename(tree1) + ".clean"
  String tree2_cleaned = basename(tree2) + ".clean"
  command <<<
    # set -euo pipefail to avoid silent failure
    set -euo pipefail

    # grab the theiaphylo version
    phylovalidate --version | tee VERSION

    # clean the trees, report if they are bifurcating
    Rscript /theiaphylo/theiaphylo/clean_phylo.R ~{tree1} > ~{tree1_cleaned} 2> >(cut -f 2 -d ' ' > TREE1_BIFURCATING)
    Rscript /theiaphylo/theiaphylo/clean_phylo.R ~{tree2} > ~{tree2_cleaned} 2> >(cut -f 2 -d ' ' > TREE2_BIFURCATING)

    # set bash variables to check them for population in conditionals
    max_distance=~{max_distance}

    # run comparison
    phylovalidate ~{tree1_cleaned} ~{tree2_cleaned} \
        ~{if (resolve_tip_discrepancies) then "-r" else ""} \
        --debug \
        2> >(tee -a PHYLOVALIDATE_STDERR >&2)

    # generate a cophylogeny plot
    Rscript /theiaphylo/theiaphylo/gen_cophylo.R ~{tree1_cleaned} ~{tree2_cleaned}

    # extract errors while maintaining a 0 exit code
    grep -Po "ERROR.*" PHYLOVALIDATE_STDERR > PHYLOVALIDATE_ERRORS || true

    # extract the distance
    tail -1 phylo_distances.txt | cut -f 2 | tr -d ' ' > PHYLOVALIDATE_DISTANCE

    # populate flag
    python3 <<CODE
    with open('TREE1_BIFURCATING', 'r') as f:
      tree1_bifurcating = f.read().strip()
    with open('TREE2_BIFURCATING', 'r') as f:
      tree2_bifurcating = f.read().strip()
    phylovalidate_flags = []
    if tree1_bifurcating == 'FALSE' or tree2_bifurcating == 'FALSE':
      phylovalidate_flags.append('polytomy')
    with open('PHYLOVALIDATE_ERRORS', 'r') as f:
      errors = set(x.strip() for x in f)
      if "ERROR - Error comparing trees: number of edges must be equal" in errors:
        phylovalidate_flags.append('edge_count_mismatch')
    with open('PHYLOVALIDATE_FLAG', 'w') as out:
      out.write(', '.join(phylovalidate_flags))
    CODE

    # run the validation
    if [[ -z ${max_distance} ]]; then
      echo "NA" > PHYLOVALIDATE
    else
      python3 <<CODE
    try:
      # check if the distance is greater than the max distance
      with open('PHYLOVALIDATE_DISTANCE', 'r') as f:
        observed_distance_str = f.read().strip()
      observed_distance = float(observed_distance_str)
      if observed_distance > ~{max_distance}:
        with open('PHYLOVALIDATE', 'w') as out:
          out.write('FAIL')
      else:
        with open('PHYLOVALIDATE', 'w') as out:
          out.write('PASS')
    # indicates that the distance is not a float, likely a None
    except ValueError:
      with open('PHYLOVALIDATE_DISTANCE', 'w') as out:
        out.write('>0')
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
    String phylovalidate_version = read_string("VERSION")
    File summary_report = "phylo_distances.txt"
    File tree1_clean = "~{tree1_cleaned}"
    File tree2_clean = "~{tree2_cleaned}"
    File cophylo_plot = "cophylo_plot.pdf"
    String phylovalidate_distance = read_string("PHYLOVALIDATE_DISTANCE")
    String phylovalidate_validation = read_string("PHYLOVALIDATE")
    String phylovalidate_flag = read_string("PHYLOVALIDATE_FLAG")
  }
}
