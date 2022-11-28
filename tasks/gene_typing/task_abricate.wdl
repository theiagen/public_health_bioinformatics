version 1.0

task abricate {
  input {
    File assembly
    String samplename
    String database
    # Parameters 
    #  --minid Minimum DNA %identity [80]
    # --mincov Minimum DNA %coverage [80]
    Int? minid
    Int? mincov
    Int cpu = 2
    String docker = "staphb/abricate:1.0.1-abaum-plasmid"
  }
  command <<<
    date | tee DATE
    abricate -v | tee ABRICATE_VERSION
    abricate --list
    abricate --check
    
    abricate \
      --db ~{database} \
      ~{'--minid ' + minid} \
      ~{'--mincov ' + mincov} \
      --threads ~{cpu} \
      --nopath \
      ~{assembly} > ~{samplename}_abricate_hits.tsv
    
    # parse out gene names into list of strings, comma-separated, final comma at end removed by sed
    abricate_genes=$(awk -F '\t' '{ print $6 }' ~{samplename}_abricate_hits.tsv | tail -n+2 | tr '\n' ',' | sed 's/.$//')

    # if variable for list of genes is EMPTY, write string saying it is empty to float to Terra table
    if [ -z "${abricate_genes}" ]; then
       abricate_genes="No genes detected by ABRicate"
    fi

    # create final output strings
    echo "${abricate_genes}" > ABRICATE_GENES
  >>>
  output {
    File abricate_results = "~{samplename}_abricate_hits.tsv"
    String abricate_genes = read_string("ABRICATE_GENES")
    String abricate_database = database
    String abricate_version = read_string("ABRICATE_VERSION")
    String abricate_docker = docker 
  }
  runtime {
    memory: "8 GB"
    cpu: cpu
    docker: docker
    disks: "local-disk 100 HDD"
  }
}
