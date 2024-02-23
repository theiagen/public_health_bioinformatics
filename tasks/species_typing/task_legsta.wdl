version 1.0

task legsta {
  meta {
    description: "Typing of Legionella pneumophila assemblies"
  }
  input {
    File assembly
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/legsta:0.5.1--hdfd78af_2"
    Int disk_size = 100
    Int memory = 8
    Int cpu = 2
  }
  command <<<
    echo $(legsta --version 2>&1) | sed 's/^.*legsta //; s/ .*\$//;' | tee VERSION
    
    legsta \
      ~{assembly} > ~{samplename}.tsv
    
    # parse outputs
    if [ ! -f ~{samplename}.tsv ]; then
      SBT="No SBT predicted"
    else
      SBT="ST$(tail -n 1 ~{samplename}.tsv | cut -f 2)"
        if [ "$SBT" == "ST-" ]; then
          SBT="No SBT predicted"
        else
          if [ "$SBT" == "ST" ]; then
            SBT="No SBT predicted"
          fi
        fi  
    fi

    echo $SBT | tee LEGSTA_SBT

  >>>
  output {
    File legsta_results = "~{samplename}.tsv"
    String legsta_predicted_sbt = read_string("LEGSTA_SBT")
    String legsta_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}
