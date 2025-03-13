version 1.0

task phylovalidate {
  input {
    File tree1_path
    File tree2_path

    Boolean? unrooted = true
    Float? rf_max_distance
    
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/theiavalidate:0.1.0"  # update!!!
    Int disk_size = 10
    Int memory = 4
    Int cpu = 1
  }
  meta {
    # added so that call caching is always turned off
    volatile: true
  }
  command <<<
    # grab the ete3 version
    ete3 --version > VERSION

    # run ete3 compare
    if ~{unrooted}; then
      ete3 compare -t ~{tree1_path} -r ~{tree2_path} --unrooted > phylocompare.txt
    else
      ete3 compare -t ~{tree1_path} -r ~{tree2_path} > phylocompare.txt
    fi

    # extract the RF distance
    tail -1 phylocompare.txt | cut -f 5 -d '|' | tr -d ' ' > PHYLOCOMPARE_RF_DISTANCE

    # run the comparison
    if [ -z ~{rf_max_distance} ]; then
      echo "NA" > PHYLOVALIDATE
    else
      python3 -c "if float(open('PHYLOCOMPARE_RF_DISTANCE', 'r').read().strip()) > ~{rf_max_distance}: open('PHYLOVALIDATE', 'w').write('FAIL'); else: open('PHYLOVALIDATE', 'w').write('PASS')"
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
    String ete3_version = read_string("VERSION")
    File summary_report = "phylocompare.txt"
    Float rf_distance = read_float("PHYLOCOMPARE_RF_DISTANCE")
    String phylovalidate = read_string("PHYLOVALIDATE")
  }
}
