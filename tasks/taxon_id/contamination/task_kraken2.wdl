version 1.0

task kraken2_theiacov {
  input {
    File read1
    File? read2
    String samplename
    File kraken2_db = "gs://theiagen-large-public-files-rp/terra/databases/kraken2/kraken2_humanGRCh38_viralRefSeq_20240828.tar.gz"
    Int cpu = 4
    Int memory = 8
    String? target_organism
    Int disk_size = 100
    String docker_image = "us-docker.pkg.dev/general-theiagen/staphb/kraken2:2.1.2-no-db"
  }
  command <<<
    # date and version control
    date | tee DATE
    kraken2 --version | head -n1 | tee VERSION
    num_reads=$(ls *fastq.gz 2> /dev/nul | wc -l)

    # Decompress the Kraken2 database
    mkdir db
    tar -C ./db/ -xzvf ~{kraken2_db} 

    if ! [ -z ~{read2} ]; then
      mode="--paired"
    fi
    echo $mode

     # determine if reads are compressed
    if [[ ~{read1} == *.gz ]]; then
      echo "Reads are compressed..."
      compressed="--gzip-compressed"
    fi
    echo $compressed

    # Run Kraken2
    kraken2 $mode $compressed \
      --threads ~{cpu} \
      --db ./db/ \
      ~{read1} ~{read2} \
      --report ~{samplename}_kraken2_report.txt \
      --output ~{samplename}.classifiedreads.txt
    
    # Compress and cleanup
    gzip ~{samplename}.classifiedreads.txt

    # capture human percentage
    percentage_human=$(grep "Homo sapiens" ~{samplename}_kraken2_report.txt | cut -f 1)
    if [ -z "$percentage_human" ] ; then percentage_human="0" ; fi
    echo $percentage_human | tee PERCENT_HUMAN

    # capture target org percentage
    if [ ! -z "~{target_organism}" ]; then
      echo "Target org designated: ~{target_organism}"
      # if target organisms is sc2, report it in a special legacy column called PERCENT_SC2
      if [[ "~{target_organism}" == "Severe acute respiratory syndrome coronavirus 2" ]]; then
        percentage_sc2=$(grep "Severe acute respiratory syndrome coronavirus 2" ~{samplename}_kraken2_report.txt  | cut -f1 )
        percent_target_organism=""
        if [ -z "$percentage_sc2" ] ; then percentage_sc2="0" ; fi
      else
        percentage_sc2="" 
        percent_target_organism=$(grep "~{target_organism}" ~{samplename}_kraken2_report.txt  | cut -f1 | head -n1 )
        if [ -z "$percent_target_organism" ] ; then percent_target_organism="0" ; fi
      fi
    else
      percent_target_organism=""
      percentage_sc2=""
    fi
    echo $percentage_sc2 | tee PERCENT_SC2
    echo $percent_target_organism | tee PERCENT_TARGET_ORGANISM

  >>>
  output {
    String date = read_string("DATE")
    String version = read_string("VERSION")
    File kraken_report = "~{samplename}_kraken2_report.txt"
    Float percent_human = read_float("PERCENT_HUMAN")
    String percent_sc2 = read_string("PERCENT_SC2")
    String percent_target_organism = read_string("PERCENT_TARGET_ORGANISM")
    String? kraken_target_organism = target_organism
    File kraken2_classified_report = "~{samplename}.classifiedreads.txt.gz" 
    String docker = docker_image
    String database = kraken2_db
  }
  runtime {
    docker: docker_image
    memory: "~{memory} GB"
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
    Int memory = 32
    Int cpu = 4
    Int disk_size = 100
  }
  command <<<
    echo $(kraken2 --version 2>&1) | sed 's/^.*Kraken version //;s/ .*$//' | tee VERSION
    date | tee DATE

    # Decompress the Kraken2 database
    mkdir db
    tar -C ./db/ -xzvf ~{kraken2_db}  

    # determine if paired-end or not
    if ! [ -z ~{read2} ]; then
      echo "Reads are paired..."
      mode="--paired"
    fi

    # determine if reads are compressed
    if [[ ~{read1} == *.gz ]]; then
      echo "Reads are compressed..."
      compressed="--gzip-compressed"
    fi

    # Run Kraken2
    echo "Running Kraken2..."
    kraken2 $mode $compressed \
        --db ./db/ \
        --threads ~{cpu} \
        --report ~{samplename}.report.txt \
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
    String kraken2_database = kraken2_db
  }
  runtime {
      docker: "~{docker}"
      memory: "~{memory} GB"
      cpu: cpu
      disks: "local-disk " + disk_size + " SSD"
      preemptible: 0
  }
}

task kraken2_parse_classified {
  input {
    File kraken2_classified_report
    File kraken2_report
    String samplename
    String? target_organism
    Int cpu = 4
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-08-28-v4"
    Int memory = 32
  }
  command <<<

    gunzip -c ~{kraken2_classified_report} > ~{samplename}.classifiedreads.txt

    python3 <<CODE 
    import pandas as pd 

    # load files into dataframe for parsing
    reads_table = pd.read_csv("~{samplename}.classifiedreads.txt", names=["classified_unclassified","read_id","taxon_id","read_len","other"], header=None, delimiter="\t")
    report_table = pd.read_csv("~{kraken2_report}", names=["percent","num_reads","num_reads_with_taxon","rank","taxon_id","name"], header=None, delimiter="\t")

    # create dataframe to store results
    results_table = pd.DataFrame(columns=["percent","num_basepairs","rank","taxon_id","name"])

    # calculate total basepairs
    total_basepairs = reads_table["read_len"].sum()

    # for each taxon_id in the report table, check if exists in the classified reads output and if so get the percent, num basepairs, rank, and name
    # write results to dataframe
    for taxon_id in report_table["taxon_id"].unique():
      if taxon_id in reads_table["taxon_id"].unique():
      
        taxon_name = report_table.loc[report_table["taxon_id"] == taxon_id, "name"].values[0].strip()
        rank = report_table.loc[report_table["taxon_id"] == taxon_id, "rank"].values[0].strip()
        taxon_basepairs = reads_table.loc[reads_table["taxon_id"] == taxon_id, "read_len"].sum()
        taxon_percent = taxon_basepairs / total_basepairs * 100
        
        results_table = pd.concat([results_table, pd.DataFrame({"percent":taxon_percent, "num_basepairs":taxon_basepairs, "rank":rank, "taxon_id":taxon_id, "name":taxon_name}, index=[0])], ignore_index=True)

    # write results to file
    results_table.to_csv("~{samplename}.report_parsed.txt", sep="\t", index=False, header=False)
    CODE

    # theiacov parsing blocks - percent human, sc2 and target organism
    # capture human percentage
    percentage_human=$(grep "Homo sapiens" ~{samplename}.report_parsed.txt | cut -f 1)
    if [ -z "$percentage_human" ] ; then percentage_human="0" ; fi
    echo $percentage_human | tee PERCENT_HUMAN

    # capture target org percentage
    if [ ! -z "~{target_organism}" ]; then
      echo "Target org designated: ~{target_organism}"
      # if target organisms is sc2, report it in a special legacy column called PERCENT_SC2
      if [[ "~{target_organism}" == "Severe acute respiratory syndrome coronavirus 2" ]]; then
        percentage_sc2=$(grep "Severe acute respiratory syndrome coronavirus 2" ~{samplename}.report_parsed.txt  | cut -f1 )
        percent_target_organism=""
        if [ -z "$percentage_sc2" ] ; then percentage_sc2="0" ; fi
      else
        percentage_sc2="" 
        percent_target_organism=$(grep "~{target_organism}" ~{samplename}.report_parsed.txt  | cut -f1 | head -n1 )
        if [ -z "$percent_target_organism" ] ; then percent_target_organism="0" ; fi
      fi
    else
      percent_target_organism=""
      percentage_sc2=""
    fi
    echo $percentage_sc2 | tee PERCENT_SC2
    echo $percent_target_organism | tee PERCENT_TARGET_ORGANISM
    
  >>>
  output {
    File kraken_report = "~{samplename}.report_parsed.txt"
    Float percent_human = read_float("PERCENT_HUMAN")
    String percent_sc2 = read_string("PERCENT_SC2")
    String percent_target_organism = read_string("PERCENT_TARGET_ORGANISM")
    String? kraken_target_organism = target_organism
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 0
  }
}