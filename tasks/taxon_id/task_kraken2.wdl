version 1.0

task kraken2_theiacov {
  input {
    File read1
    File? read2
    String samplename
    String kraken2_db = "/kraken2-db"
    Int cpu = 4
    String? target_organism
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

    percentage_human=$(grep "Homo sapiens" ~{samplename}_kraken2_report.txt | cut -f 1)
     # | tee PERCENT_HUMAN
    percentage_sc2=$(grep "Severe acute respiratory syndrome coronavirus 2" ~{samplename}_kraken2_report.txt | cut -f1 )
     # | tee PERCENT_COV
    if [ -z "$percentage_human" ] ; then percentage_human="0" ; fi
    if [ -z "$percentage_sc2" ] ; then percentage_sc2="0" ; fi
    echo $percentage_human | tee PERCENT_HUMAN
    echo $percentage_sc2 | tee PERCENT_SC2
    # capture target org percentage 
    if [ ! -z "~{target_organism}" ]; then
      echo "Target org designated: ~{target_organism}"
      percent_target_organism=$(grep "~{target_organism}" ~{samplename}_kraken2_report.txt | cut -f1 | head -n1 )
      if [-z "$percent_target_organism" ] ; then percent_target_organism="0" ; fi
    else 
      percent_target_organism=""
    fi
    echo $percent_target_organism | tee PERCENT_TARGET_ORGANISM

  >>>
  output {
    String date = read_string("DATE")
    String version = read_string("VERSION")
    File kraken_report = "~{samplename}_kraken2_report.txt"
    Float percent_human = read_float("PERCENT_HUMAN")
    Float percent_sc2 = read_float("PERCENT_SC2")
    String percent_target_organism = read_string("PERCENT_TARGET_ORGANISM")
    String? kraken_target_organism = target_organism
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/kraken2:2.0.8-beta_hv"
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