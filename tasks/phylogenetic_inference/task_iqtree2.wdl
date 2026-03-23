version 1.0

task iqtree2 {
  input {
    File alignment
    String cluster_name
    String? iqtree2_model # For comparison to other tools use HKY for bactopia, GTR+F+I for grandeur, GTR+G4 for nullarbor, GTR+G for dryad
    Int iqtree2_ultrafast_bootstraps = 1000 #  Ultrafast bootstrap replicates
    Int alrt = 1000 # SH-like approximate likelihood ratio test (SH-aLRT) replicates
    String? iqtree2_opts

    String docker = "us-docker.pkg.dev/general-theiagen/staphb/iqtree2:2.1.2"
    Int disk_size = 100
    Int cpu = 4
    Int memory = 32
  }
  command <<<
    # date and version control
    date | tee DATE

    # multiple sed statements to get down to a string that is just "version 2.1.2"
    iqtree2 --version | grep version | sed 's|.*version|version|;s| COVID-edition for Linux.*||' | tee VERSION

    # check if iqtree2_model input is set and output for sanity
    if [ -n "~{iqtree2_model}" ]; then
      echo "DEBUG: User provided iqtree2_model ~{iqtree2_model}, will use this for running iqtree2"
      IQTREE2_MODEL="~{iqtree2_model}"
    else
      echo "DEBUG: User did not supply an iqtree2_model input, will use iqtree2's model finder"
    fi

    # sanity check
    echo "DEBUG: IQTREE2_MODEL is set to: " ${IQTREE2_MODEL}

    # make sure there are more than 3 genomes in the dataset
    numGenomes=$(grep -o '>' ~{alignment} | wc -l)    
    if [ "$numGenomes" -gt 3 ]; then
      cp ~{alignment} ./msa.fasta

      # run iqtree2
      #   -nt : number of CPU cores for multicore version
      #   -s : input alignment file
      #   -m : model
      #   -bb : number of bootstrap replicates
      #   -alrt : number of replicates to perform SH-like approximate likelihood ration test  
      if [[ -v IQTREE2_MODEL ]] ; then # iqtree2 model set; use -m tag
        echo "DEBUG: running iqtree2 with the -m flag which is used to provide a model; user-specified " ${IQTREE2_MODEL}
        iqtree2 \
          -nt AUTO \
          -s msa.fasta \
          -m ${IQTREE2_MODEL} \
          -bb ~{iqtree2_ultrafast_bootstraps} \
          -alrt ~{alrt} ~{iqtree2_opts}

        # write the iqtree2_model used to a txt file for output as a string
        echo ${IQTREE2_MODEL} | tee IQTREE2_MODEL.TXT

      else # iqtree model is not set; do not use -m tag
        echo "DEBUG: running iqtree2 without the -m flag which is used to provide a model. Will default to iqtree2 default (Model Finder)"
        iqtree2 \
          -nt AUTO \
          -s msa.fasta \
          -bb ~{iqtree2_ultrafast_bootstraps} \
          -alrt ~{alrt} ~{iqtree2_opts}

        # for scenario where user did not specify iqtree2_model input nor core_genome boolean input, determine iqtree2_model used by parsing log file
        # first sed is to remove "Best-fit model: " and second sed is to remove anything after the word "chosen *", leaving only the name of the model
        grep "Best-fit model" msa.fasta.log | sed 's|Best-fit model: ||g;s|chosen.*||' | tee IQTREE2_MODEL.TXT

      fi

      # rename the final output newick file
      cp -v msa.fasta.contree ~{cluster_name}_iqtree.nwk
    else
      echo "DEBUG: not enough genomes provided; more than 3 are required to run iqtree2"
    fi
  >>>
  output {
    String date = read_string("DATE")
    String iqtree2_version = read_string("VERSION")
    File ml_tree = "~{cluster_name}_iqtree.nwk"
    String iqtree2_model_used = read_string("IQTREE2_MODEL.TXT")
    String iqtree2_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}
