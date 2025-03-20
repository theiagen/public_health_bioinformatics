version 1.0

task phylovalidate {
  input {
    File tree1_path
    File tree2_path

    String? root_tips
    Boolean? unrooted
    Float? lrm_max_dist
    
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

    # run the comparison
    if [[ -z ~{root_tips} ]]; then
      phylocompare.py ~{tree1_path} ~{tree2_path} \
        ~{true="--unrooted" false="" unrooted} \
        --debug > phylocompare.txt
    else
      phylocompare.py ~{tree1_path} ~{tree2_path} \
        --root-tips ~{root_tips} \
        --debug > phylocompare.txt
    fi

    # extract the RF distance
    tail -1 phylocompare.txt | cut -f 3 | tr -d ' ' > PHYLOCOMPARE_RF_DISTANCE

    # run the comparison
    if [ -z ~{lrm_max_dist} ]; then
      echo "NA" > PHYLOVALIDATE
    else
      python3 -c "if float(open('PHYLOCOMPARE_LRM_DISTANCE', 'r').read().strip()) > ~{lrm_max_dist}: open('PHYLOVALIDATE', 'w').write('FAIL'); else: open('PHYLOVALIDATE', 'w').write('PASS')"
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
    File summary_report = "phylocompare.txt"
    Float lrm_distance = read_float("PHYLOCOMPARE_LRM_DISTANCE")
    String phylovalidate = read_string("PHYLOVALIDATE")
  }
}
