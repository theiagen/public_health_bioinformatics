version 1.0

task plasmidfinder {
  input {
    File assembly
    String samplename
    Int cpu = 8
    Int memory = 16
    String docker = "staphb/plasmidfinder:2.1.6"
    String? database
    String? database_path
    String? method_path
    # minimum coverage threshold
    Float? min_cov 
    # minimum blast identity threshold
    Float? threshold

  }
  command <<<  
  date | tee DATE

  if [[ ! -z "~{database}" ]]; then 
    echo "User database identified; ~{database} will be utilized for analysis"
    plasmidfinder_db_version="~{database}"
  else
    plasmidfinder_db_version="unmodified from plasmidfinder docker container"
  fi

  echo ${plasmidfinder_db_version} | tee PLASMIDFINDER_DB_VERSION

  plasmidfinder.py \
  -i ~{assembly} \
  -x \
  ~{'-d ' + database} \
  ~{'-p ' + database_path} \
  ~{'-mp ' + method_path} \
  ~{'-l ' + min_cov} \
  ~{'-t ' + threshold} 

  # parse outputs
  if [ ! -f results_tab.tsv ]; then
    PF="No plasmids detected in database"
  else
    PF="$(tail -n +2 results_tab.tsv | cut -f 2 | sort | uniq -u | paste -s -d, - )"
      if [ "$PF" == "" ]; then
        PF="No plasmids detected in database"
      fi  
  fi
  echo $PF | tee PLASMIDS

  mv results_tab.tsv ~{samplename}_results.tsv
  mv Hit_in_genome_seq.fsa ~{samplename}_seqs.fsa

  >>>
  output {
    String plasmidfinder_plasmids = read_string("PLASMIDS")
    File plasmidfinder_results = "~{samplename}_results.tsv"
    File plasmidfinder_seqs = "~{samplename}_seqs.fsa"
    String plasmidfinder_docker = docker
    String plasmidfinder_db_version = read_string("PLASMIDFINDER_DB_VERSION")
  }
  runtime {
    memory: "~{memory} GB"
    cpu: cpu
    docker: "~{docker}"
    disks: "local-disk 100 HDD"
  }
}
