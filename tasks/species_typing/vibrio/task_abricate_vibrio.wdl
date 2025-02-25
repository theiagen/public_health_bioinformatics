version 1.0

task abricate_vibrio {
  input {
    File assembly
    String samplename
    String database = "vibrio"
    Int min_percent_identity
    Int min_percent_coverage
    Int cpu = 2
    Int memory = 4
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/abricate:1.0.1-vibrio-cholera"
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
    
    # presence or absence genes - ctxA, ompW and toxR
    # if empty report as (not detected)
    grep -E 'ctxA' ~{samplename}_abricate_hits.tsv > ctxA
    if [ ! -s ctxA ]; then echo "(not detected)" > ctxA; else echo "present" > ctxA; fi
    grep -E 'ompW' ~{samplename}_abricate_hits.tsv > ompW
    if [ ! -s ompW ]; then echo "(not detected)" > ompW; else echo "present" > ompW; fi
    grep -E 'toxR' ~{samplename}_abricate_hits.tsv > toxR
    if [ ! -s toxR ]; then echo "(not detected)" > toxR; else echo "present" > toxR; fi

    # biotype - tcpA classical or tcpA ElTor
    grep -E 'tcpA_classical' ~{samplename}_abricate_hits.tsv > tcpA_classical
    grep -E 'tcpA_ElTor' ~{samplename}_abricate_hits.tsv > tcpA_ElTor
    if [ ! -s tcpA_classical ] && [ ! -s tcpA_ElTor ]; then 
      echo "(not detected)" > BIOTYPE; 
    else 
      if [ ! -s tcpA_classical ] && [ -s tcpA_ElTor ]; then 
        echo "tcpA_ElTor" > BIOTYPE; 
      else 
        echo 'tcpA_Classical' > BIOTYPE
      fi
    fi

    # serogroup - O1 or O139
    grep -E 'wbeN_O1' ~{samplename}_abricate_hits.tsv > wbeN_O1
    grep -E 'wbfR_O139' ~{samplename}_abricate_hits.tsv > wbfR_O139
    if [ ! -s wbeN_O1 ] && [ ! -s wbfR_O139 ]; then 
      echo "(not detected)" > SEROGROUP; 
    else 
      if [ ! -s wbeN_O1 ] && [ -s wbfR_O139 ]; then 
        echo "O139" > SEROGROUP; 
      else 
        echo 'O1' > SEROGROUP
      fi
    fi
  >>>
  output {
    File abricate_vibrio_results = "~{samplename}_abricate_hits.tsv"
    String abricate_vibrio_database = database
    String abricate_vibrio_docker = docker
    String abricate_vibrio_version = read_string("ABRICATE_VERSION")
    String abricate_vibrio_ctxA = read_string("ctxA")
    String abricate_vibrio_ompW = read_string("ompW")
    String abricate_vibrio_toxR = read_string("toxR")
    String abricate_vibrio_biotype = read_string("BIOTYPE")
    String abricate_vibrio_serogroup = read_string("SEROGROUP")
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
