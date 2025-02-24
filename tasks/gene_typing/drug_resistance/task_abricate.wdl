version 1.0

task abricate {
  input {
    File assembly
    String samplename
    String database
    # Parameters 
    #  --minid Minimum DNA %identity [80]
    # --mincov Minimum DNA %coverage [80]
    Int? min_percent_identity
    Int? min_percent_coverage
    Int cpu = 2
    Int memory = 8
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/abricate:1.0.1-abaum-plasmid"
  }
  command <<<
    date | tee DATE
    abricate -v | tee ABRICATE_VERSION
    abricate --list
    abricate --check
    
    abricate \
      --db ~{database} \
      ~{'--minid ' + min_percent_identity} \
      ~{'--mincov ' + min_percent_coverage} \
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
    memory: memory + " GB"
    cpu: cpu
    docker: docker
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
  }
}

task abricate_flu {
  input {
    File assembly
    String samplename
    String database = "insaflu"
    Int min_percent_identity = 70
    Int min_percent_coverage = 60
    Int cpu = 2
    Int memory = 4
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/abricate:1.0.1-insaflu-220727"
    Int disk_size = 100
  }
  command <<<
    date | tee DATE    
    abricate -v | tee ABRICATE_VERSION
    # run abricate
    abricate \
      --db ~{database} \
      ~{'--minid ' + min_percent_identity} \
      ~{'--mincov ' + min_percent_coverage} \
      --threads ~{cpu} \
      --nopath \
      ~{assembly} > ~{samplename}_abricate_hits.tsv

    # capturing flu type (A or B based on M1 hit) and subtype (e.g. H1 and N1 based on HA/NA hits)
    ## awk for gene column ($6) to grab subtype ($15)
    cat ~{samplename}_abricate_hits.tsv | awk -F '\t' '{if ($6=="M1") print $15}' > FLU_TYPE
    HA_hit=$(cat ~{samplename}_abricate_hits.tsv | awk -F '\t' '{if ($6=="HA") print $15 }')
    NA_hit=$(cat ~{samplename}_abricate_hits.tsv | awk -F '\t' '{if ($6=="NA") print $15 }')
    if [[ ! (-z "${HA_hit}")  &&  ! (-z "${NA_hit}") ]]; then
      flu_subtype="${HA_hit}${NA_hit}" && echo "$flu_subtype" >  FLU_SUBTYPE
    fi
    if [[ -z "${HA_hit}" ]]; then
      flu_subtype="${NA_hit}" && echo "$flu_subtype" >  FLU_SUBTYPE
    elif [[ -z "${NA_hit}" ]]; then
      flu_subtype="${HA_hit}" && echo "$flu_subtype" >  FLU_SUBTYPE
    else
      flu_subtype="${HA_hit}${NA_hit}" && echo "$flu_subtype" >  FLU_SUBTYPE
    fi
    #flu_subtype="${HA_hit}${NA_hit}" && echo "$flu_subtype" >  FLU_SUBTYPE
  >>>
  output {
    String abricate_flu_type = read_string("FLU_TYPE")
    String abricate_flu_subtype = read_string("FLU_SUBTYPE")
    File abricate_flu_results = "~{samplename}_abricate_hits.tsv"
    String abricate_flu_database = database
    String abricate_flu_version = read_string("ABRICATE_VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible:  0
  }
}