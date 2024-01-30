version 1.0

task kraken2_theiacov {
  input {
    File read1
    File? read2
    String samplename
    String? organism
    String kraken2_db = "/kraken2-db"
    Int cpu = 4
    String? target_org
    Int disk_size = 100
  }
  command <<<
    # date and version control
    date | tee DATE
    kraken2 --version | head -n1 | tee VERSION
    num_reads=$(ls *fastq.gz 2> /dev/nul | wc -l)
    if ! [ -z ~{read2} ]; then
      mode="--paired"
    fi
    echo $mode
    kraken2 $mode \
      --threads ~{cpu} \
      --db ~{kraken2_db} \
      ~{read1} ~{read2} \
      --report ~{samplename}_kraken2_report.txt >/dev/null

    # capture percentage human
    percentage_human=$(grep "Homo sapiens" ~{samplename}_kraken2_report.txt | cut -f 1)
    if [ -z "$percentage_human" ] ; then percentage_human="0" ; fi
    echo $percentage_human | tee PERCENT_HUMAN

    # Capture percentage of SARS-CoV-2, regardless of organism input
    percentage_sc2=$(grep "Severe acute respiratory syndrome coronavirus 2" ~{samplename}_kraken2_report.txt | cut -f1 )
    if [ -z "$percentage_sc2" ] ; then percentage_sc2="0" ; fi
    echo $percentage_sc2 | tee PERCENT_SC2

    # capture percentage of various TheiaCoV pathogens based on the organism input
    if [[ "~{organism}" == "sars-cov-2" ]] ; then
      kraken_output=$(grep "Severe acute respiratory syndrome coronavirus" ~{samplename}_kraken2_report.txt | sort -k2,2n | tail -1)
      kraken_identified_organism=$(echo "$kraken_output" | awk '{ $1=$2=$3=$4=$5=""; print $0 }')
      kraken_identified_organism_abundance=$(echo "$kraken_output" | cut -f1)
      # for legacy purposes
      echo $kraken_identified_organism_abundance > PERCENT_SC2
    elif [[ "~{organism}" == "flu" ]] ; then
      kraken_output=$(grep "Influenza" ~{samplename}_kraken2_report.txt | sort -k3,3n | tail -1)
      kraken_identified_organism=$(echo "$kraken_output" | awk '{ $1=$2=$3=$4=$5=""; print $0 }')
      kraken_identified_organism_abundance=$(echo "$kraken_output" | cut -f1)
    elif [[ "~{organism}" == "MPXV" ]] ; then
      kraken_output=$(grep "Monkeypox" ~{samplename}_kraken2_report.txt | sort -k3,3n | tail -1)
      kraken_identified_organism=$(echo "$kraken_output" | awk '{ $1=$2=$3=$4=$5=""; print $0 }')
      kraken_identified_organism_abundance=$(echo "$kraken_output" | cut -f1)
    elif [[ "~{organism}" == "HIV" ]] ; then
      kraken_output=$(grep "Human immunodeficiency virus" ~{samplename}_kraken2_report.txt | sort -k3,3n | tail -1)
      kraken_identified_organism=$(echo "$kraken_output" | awk '{ $1=$2=$3=$4=$5=""; print $0 }')
      kraken_identified_organism_abundance=$(echo "$kraken_output" | cut -f1)
    elif [[ "~{organism}" == "WNV" ]] ; then
      kraken_output=$(grep "West Nile virus" ~{samplename}_kraken2_report.txt | sort -k3,3n | tail -1)
      kraken_identified_organism=$(echo "$kraken_output" | awk '{ $1=$2=$3=$4=$5=""; print $0 }')
      kraken_identified_organism_abundance=$(echo "$kraken_output" | cut -f1)
    else
      kraken_identified_organism=""
      kraken_identified_organism_abundance=""
    fi

    echo $kraken_identified_organism | tee MOST_ABUNDANT_ORGANISM
    echo $kraken_identified_organism_abundance | tee PERCENT_MOST_ABUNDANT_ORGANISM

    # capture target org percentage 
    if [ ! -z "~{target_org}" ]; then 
      echo "Target org designated: ~{target_org}"
      percent_target_org=$(grep "~{target_org}" ~{samplename}_kraken2_report.txt | cut -f1 | head -n1 )
      if [-z "$percent_target_org" ] ; then percent_target_org="0" ; fi
    else 
      percent_target_org=""
    fi
    echo $percent_target_org | tee PERCENT_TARGET_ORG
  >>>
  output {
    String date = read_string("DATE")
    String version = read_string("VERSION")
    File kraken_report = "~{samplename}_kraken2_report.txt"
    Float percent_human = read_float("PERCENT_HUMAN")
    Float percent_sc2 = read_float("PERCENT_SC2")
    String percent_target_org = read_string("PERCENT_TARGET_ORG")
    String? kraken_target_org = target_org
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/kraken2:2.0.8-beta_hv"
    String most_abundant_organism = read_string("MOST_ABUNDANT_ORGANISM")
    String percent_most_abundant_organism = read_string("PERCENT_MOST_ABUNDANT_ORGANISM")
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/staphb/kraken2:2.0.8-beta_hv"
    memory: "8 GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}

# standalone version (no default database included)
task kraken2_standalone {
  input {
    File read1
    File? read2
    File kraken2_db
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/kraken2:2.1.2-no-db"
    String kraken2_args = ""
    String classified_out = "classified#.fastq"
    String unclassified_out = "unclassified#.fastq"
    Int mem = 32
    Int cpu = 4
  }
  command <<<
    echo $(kraken2 --version 2>&1) | sed 's/^.*Kraken version //;s/ .*$//' | tee VERSION
    date | tee DATE

    # Decompress the Kraken2 database
    mkdir db
    tar -C ./db/ -xzvf ~{kraken2_db}  

    # determine if paired-end or not
    if ! [ -z ~{read2} ]; then
      mode="--paired"
    fi

    # Run Kraken2
    kraken2 $mode \
        --db ./db/ \
        --threads ~{cpu} \
        --report ~{samplename}.report.txt \
        --gzip-compressed \
        --unclassified-out ~{samplename}.~{unclassified_out} \
        --classified-out ~{samplename}.~{classified_out} \
        --output ~{samplename}.classifiedreads.txt \
        ~{kraken2_args} \
        ~{read1} ~{read2}
    
    # Compress and cleanup
    gzip *.fastq
    gzip ~{samplename}.classifiedreads.txt

    # Report percentage of human reads
    percentage_human=$(grep "Homo sapiens" ~{samplename}.report.txt | cut -f 1)
    if [ -z "$percentage_human" ] ; then percentage_human="0" ; fi
    echo $percentage_human | tee PERCENT_HUMAN
    
    # rename classified and unclassified read files if SE
    if [ -e "~{samplename}.classified#.fastq.gz" ]; then
      mv "~{samplename}.classified#.fastq.gz" ~{samplename}.classified_1.fastq.gz
    fi
    if [ -e "~{samplename}.unclassified#.fastq.gz" ]; then
      mv "~{samplename}.unclassified#.fastq.gz" ~{samplename}.unclassified_1.fastq.gz
    fi

  >>>
  output {
    String kraken2_version = read_string("VERSION")
    String kraken2_docker = docker
    String analysis_date = read_string("DATE")
    File kraken2_report = "~{samplename}.report.txt"
    File kraken2_classified_report = "~{samplename}.classifiedreads.txt.gz"
    File kraken2_unclassified_read1 = "~{samplename}.unclassified_1.fastq.gz"
    File? kraken2_unclassified_read2 = "~{samplename}.unclassified_2.fastq.gz"
    File kraken2_classified_read1 = "~{samplename}.classified_1.fastq.gz"
    Float kraken2_percent_human = read_float("PERCENT_HUMAN")
    File? kraken2_classified_read2 = "~{samplename}.classified_2.fastq.gz"
  }
  runtime {
      docker: "~{docker}"
      memory: "~{mem} GB"
      cpu: cpu
      disks: "local-disk 100 SSD"
      preemptible: 0
  }
}