version 1.0

task iqtree2 {
  input {
    File alignment
    String cluster_name
    String? iqtree2_model # For comparison to other tools use HKY for bactopia, GTR+F+I for grandeur, GTR+G4 for nullarbor, GTR+G for dryad
    Int iqtree2_bootstraps = 1000 #  Ultrafast bootstrap replicates
    Int alrt = 1000 # SH-like approximate likelihood ratio test (SH-aLRT) replicates
    String? iqtree2_opts
    String docker = "quay.io/staphb/iqtree2:2.1.2"
    Int disk_size = 100
    Int cpu = 4
    Int memory = 32
    Boolean? core_genome
  }
  command <<<
    # date and version control
    date | tee DATE
    # multiple sed statements to get down to a string that is just "version 2.1.2"
    iqtree2 --version | grep version | sed 's|.*version|version|;s| COVID-edition for Linux.*||' | tee VERSION

    # if iqtree2_model is set by user, use that String input
    if [ -n "~{iqtree2_model}" ]; then
      echo "user provided an iqtree2_model string input, will use this for running iqtree2"
      IQTREE2_MODEL="~{iqtree2_model}"
    else
      echo "User did not supply an iqtree2_model input, setting based on boolean core_genome"

      # if iqtree2_model is NOT set by user, then set iqtree2_model based on boolean core_genome
      # if core_genome is set to TRUE, then use model "GTR+G"
      if [[ "~{core_genome}" == true ]]; then
        echo "core_genome boolean was set to true, so using iqtree2_model GTR+G"
        IQTREE2_MODEL="GTR+G"
      elif [ "~{core_genome}" == false ]; then
        echo "core_genome boolean was set to false, so using iqtree2_model GTR+I+G"
        IQTREE2_MODEL="GTR+I+G"
      else
        echo "iqtree2_model was not specified by user AND core_genome was not specified, so we will use the default setting from iqtree2"
      fi
    fi

    # sanity check
    echo "IQTREE2_MODEL is set to:" ${IQTREE2_MODEL}

    numGenomes=`grep -o '>' ~{alignment} | wc -l`
    if [ "$numGenomes" -gt 3 ]
    then
      cp ~{alignment} ./msa.fasta

      # run iqtree2
      # if IQTREE2_MODEL bash variable is set, use -m flag, otherwise do not use -m flag
      if [[ -v IQTREE2_MODEL ]] ; then
        iqtree2 \
          -nt AUTO \
          -s msa.fasta \
          -m ${IQTREE2_MODEL} \
          -bb ~{iqtree2_bootstraps} \
          -alrt ~{alrt} ~{iqtree2_opts}
      else
        echo "running iqtree2 without the -m flag for providing a model. Will default to iqtree2 defaults"
        iqtree2 \
          -nt AUTO \
          -s msa.fasta \
          -bb ~{iqtree2_bootstraps} \
          -alrt ~{alrt} ~{iqtree2_opts}
      fi

      # rename the final output newick file
      cp -v msa.fasta.contree ~{cluster_name}_iqtree.nwk
    fi
  >>>
  output {
    String date = read_string("DATE")
    String version = read_string("VERSION")
    File ml_tree = "~{cluster_name}_iqtree.nwk"
    String iqtree2_docker = docker
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
