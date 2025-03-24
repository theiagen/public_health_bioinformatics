version 1.0

task plasmidfinder {
  input {
    File assembly
    String samplename
    Int cpu = 2
    Int memory = 8
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/plasmidfinder:2.1.6"
    Int disk_size = 50
    String? database
    String? database_path
    String? method_path
    # minimum coverage threshold
    Float? min_percent_coverage 
    # minimum blast identity threshold
    Float? min_percent_identity
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
  ~{'-l ' + min_percent_coverage} \
  ~{'-t ' + min_percent_identity} 

  # parse outputs
  if [ ! -f results_tab.tsv ]; then
    PF="No plasmids detected in database"
  else
    PF="$(tail -n +2 results_tab.tsv | uniq | cut -f 2 | sort | paste -s -d, - )"
      if [ "$PF" == "" ]; then
        PF="No plasmids detected in database"
      fi  
  fi
  echo "$PF" | tee PLASMIDS

  mv -v results_tab.tsv ~{samplename}_results.tsv
  mv -v Hit_in_genome_seq.fsa ~{samplename}_seqs.fsa

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
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
  }
}
