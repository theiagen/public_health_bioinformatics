version 1.0

task iqtree2 {
  input {
    File alignment
    String cluster_name
    String iqtree_model = "GTR+G" # For comparison to other tools use HKY for bactopia, GTR+F+I for grandeur, GTR+G4 for nullarbor, GTR+G for dryad
    Int iqtree_bootstraps = 1000 #  Ultrafast bootstrap replicates
    Int alrt = 1000 # SH-like approximate likelihood ratio test (SH-aLRT) replicates
    String iqtree_opts = ""
    String docker = "quay.io/staphb/iqtree2:2.1.2"
    Int disk_size = 100
    Int cpu = 4
    Int memory = 32
  }
  command <<<
    # date and version control
    date | tee DATE
    # multiple sed statements to get down to a string that is just "version 2.1.2"
    iqtree2 --version | grep version | sed 's|.*version|version|;s| COVID-edition for Linux.*||' | tee VERSION

    numGenomes=`grep -o '>' ~{alignment} | wc -l`
    if [ "$numGenomes" -gt 3 ]
    then
      cp ~{alignment} ./msa.fasta
      iqtree2 \
      -nt AUTO \
      -s msa.fasta \
      -m ~{iqtree_model} \
      -bb ~{iqtree_bootstraps} \
      -alrt ~{alrt} \
      ~{iqtree_opts}

      cp -v msa.fasta.contree ~{cluster_name}_iqtree.nwk
    fi
  >>>
  output {
    String date = read_string("DATE")
    String version = read_string("VERSION")
    File ml_tree = "~{cluster_name}_iqtree.nwk"
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}
