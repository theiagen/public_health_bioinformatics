version 1.0

task kraken2 {
  input {
    File read1
    File? read2
    String samplename
    File kraken2_db = "gs://theiagen-public-resources-rp/reference_data/databases/kraken2/k2_viral-refseq_human-GRCh38_20260220.tar.gz"
    String kraken2_args = ""
    String classified_out = "classified#.fastq"
    String unclassified_out = "unclassified#.fastq"
    String? target_organism
    Boolean call_bracken = true
    Int? bracken_kmer_length
    Int cpu = 4
    Int memory = 32
    Int disk_size = 100
    String docker_image = "us-docker.pkg.dev/general-theiagen/theiagen/kraken2:2.17.1"
  }
  command <<<
    # fail hard
    set -eu

    # date and version control
    echo "INFO: Kraken2 version:"
    echo $(kraken2 --version 2>&1) | sed 's/^.*Kraken version //;s/ .*$//' | tee KRAKEN2_VERSION
    num_reads=$(ls *fastq.gz 2> /dev/null | wc -l)

    # Decompress the Kraken2 database
    mkdir db
    tar -C ./db/ -xzvf ~{kraken2_db} 

    if ! [ -z ~{read2} ]; then
      echo "DEBUG: Reads are paired..."
      mode="--paired"
    else
      echo "DEBUG: Reads are single-end..."
      mode=""
    fi

     # determine if reads are compressed
    if [[ ~{read1} == *.gz ]]; then
      echo "DEBUG: Reads are compressed..."
      compressed="--gzip-compressed"
    else
      echo "DEBUG: Reads are not compressed..."
      compressed=""
    fi

    # Run Kraken2
    echo "INFO: Running Kraken2..."
    kraken2 $mode $compressed \
      --threads ~{cpu} \
      --db ./db/ \
      --unclassified-out ~{samplename}.~{unclassified_out} \
      --classified-out ~{samplename}.~{classified_out} \
      --report ~{samplename}_kraken2_report.txt \
      --output ~{samplename}.classifiedreads.txt \
      ~{kraken2_args} \
      ~{read1} ~{read2}

    # Compress and cleanup
    gzip *.fastq
    gzip ~{samplename}.classifiedreads.txt

    # Run Bracken 
    touch BRACKEN_VERSION 
    if [ "~{call_bracken}" == "true" ]; then
      # check if kraken database is compatible with bracken (i.e. has kmer distribution files)
      if [ -z "$(ls db/database*mers\.kmer_distrib 2> /dev/null)" ]; then
        echo "ERROR: Bracken kmer distribution files not found in the Kraken2 database. Skipping" >&2
      else
        # if bracken_kmer_length isn't provided, infer as the kmer length directly under the mean read length
        if [ -z "~{bracken_kmer_length}" ]; then

          echo "INFO: Bracken version:"
          bracken -v | sed 's/^Bracken //' | tee BRACKEN_VERSION

          seqkit stats ~{read1} > read_lengths.tsv
          python3 <<CODE
    import os
    import re
    with open("read_lengths.tsv", "r") as f:
      # last line has the read data
      for line in f:
        data = line.strip().split()
    mean_len = round(float(data[6].replace(",",""))) # handle comma in mean read length for large numbers
    print(f"DEBUG: Inferred mean read length: {mean_len}")
    # obtain the kmer lengths available
    kmer_files = [f for f in os.listdir("db/") if f.endswith(".kmer_distrib")]
    kmer_dists = [int(re.match(r"database(\d+)mers\.kmer_distrib", f).group(1)) for f in kmer_files]
    # the best kmer is directly under the mean length size
    best_kmer = max([k for k in kmer_dists if k < mean_len], default=min(kmer_dists))
    print(f"INFO: Selected k-mer length for Bracken: {best_kmer}")
    with open("bracken_kmer_length", "w") as out:
      out.write(str(best_kmer))
    CODE

          bracken_kmer_length=$(cat bracken_kmer_length)
        else
          bracken_kmer_length="~{bracken_kmer_length}"
          echo "INFO: Using provided Bracken read length: ${bracken_kmer_length}"
        fi

        bracken -d ./db/ \
          -i ~{samplename}_kraken2_report.txt \
          -o ~{samplename}_bracken_summary.txt \
          -w ~{samplename}_bracken_report.txt \
          -r ${bracken_kmer_length} \
          -l S || true
      fi
    fi

    # Report percentage of human reads
    if [ -f "~{samplename}_bracken_report.txt" ]; then
      kraken2_report="~{samplename}_bracken_report.txt"
    else
      kraken2_report="~{samplename}_kraken2_report.txt"
    fi
    percentage_human=$(grep "Homo sapiens" $kraken2_report | cut -f 1)
    if [ -z "$percentage_human" ] ; then percentage_human="0" ; fi
    echo "INFO: Percentage human:"
    echo $percentage_human | tee PERCENT_HUMAN

    # rename classified and unclassified read files if SE
    if [ -e "~{samplename}.classified#.fastq.gz" ]; then
      mv "~{samplename}.classified#.fastq.gz" ~{samplename}.classified_1.fastq.gz
    fi
    if [ -e "~{samplename}.unclassified#.fastq.gz" ]; then
      mv "~{samplename}.unclassified#.fastq.gz" ~{samplename}.unclassified_1.fastq.gz
    fi

    # capture target org percentage
    if [ ! -z "~{target_organism}" ]; then
      echo "DEBUG: Target org designated: ~{target_organism}"
      # if target organisms is sc2, report it in a special legacy column called PERCENT_SC2
      if [[ "~{target_organism}" == "Severe acute respiratory syndrome coronavirus 2" ]]; then
        percentage_sc2=$(grep "Severe acute respiratory syndrome coronavirus 2" $kraken2_report  | cut -f1 )
        percent_target_organism=""
        if [ -z "$percentage_sc2" ]; then 
          percentage_sc2="0"
        fi
        echo "INFO: Percentage SARS-CoV-2:"
        echo "" > PERCENT_TARGET_ORGANISM
      else
        echo "" > PERCENT_SC2 
        percent_target_organism=$(grep "~{target_organism}" $kraken2_report  | cut -f1 | head -n1 )
        if [ -z "$percent_target_organism" ]; then 
          percent_target_organism="0"
        fi
        echo "INFO: Percentage target organism (~{target_organism}):"
        echo $percent_target_organism | tee PERCENT_TARGET_ORGANISM
      fi
    else
      echo "" > PERCENT_SC2
      echo "" > PERCENT_TARGET_ORGANISM
    fi
  >>>
  output {
    String kraken2_version = read_string("KRAKEN2_VERSION")
    String bracken_version = read_string("BRACKEN_VERSION")
    File kraken2_report = "~{samplename}_kraken2_report.txt"
    File? bracken_report = "~{samplename}_bracken_report.txt"
    Float kraken2_percent_human = read_float("PERCENT_HUMAN")
    String kraken2_percent_sc2 = read_string("PERCENT_SC2")
    String kraken2_percent_target_organism = read_string("PERCENT_TARGET_ORGANISM")
    String? kraken2_target_organism = target_organism
    File kraken2_classified_report = "~{samplename}.classifiedreads.txt.gz"
    File kraken2_unclassified_read1 = "~{samplename}.unclassified_1.fastq.gz"
    File? kraken2_unclassified_read2 = "~{samplename}.unclassified_2.fastq.gz"
    File kraken2_classified_read1 = "~{samplename}.classified_1.fastq.gz"
    File? kraken2_classified_read2 = "~{samplename}.classified_2.fastq.gz"
    String kraken2_database = kraken2_db 
    String docker = docker_image
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

task kraken2_parse_classified {
  input {
    File kraken2_classified_report
    File kraken2_report
    String samplename
    String? target_organism
    Int cpu = 4
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-08-28-v4"
    Int memory = 8
  }
  command <<<

    gunzip -c ~{kraken2_classified_report} > ~{samplename}.classifiedreads.txt

    python3 <<CODE 
    import csv
    from collections import defaultdict

    # Define file paths using the template variables
    reads_file = "~{samplename}.classifiedreads.txt"
    report_file = "~{kraken2_report}"
    output_file = "~{samplename}.report_parsed.txt"

    # First pass: read report file and build taxon info dictionary
    taxon_info = {}
    with open(report_file, 'r') as f:
      for line in csv.reader(f, delimiter='\t'):
        if len(line) >= 6:
          percent, num_reads, num_reads_with_taxon, rank, taxon_id, name = line
          taxon_info[taxon_id] = {'rank': rank.strip(), 'name': name.strip()}

    # Second pass: process reads file and count basepairs
    taxon_basepairs = defaultdict(int)
    total_basepairs = 0
    present_taxa = set()

    with open(reads_file, 'r') as f:
      for line in csv.reader(f, delimiter='\t'):
        if len(line) >= 5:
          _, _, taxon_id, read_len, _ = line[0], line[1], line[2], line[3], line[4]
          try:
            read_len = int(read_len)
            total_basepairs += read_len
            taxon_basepairs[taxon_id] += read_len
            present_taxa.add(taxon_id)
          except ValueError:
            continue

    # Write results to file - following the original output format
    with open(output_file, 'w', newline='') as f:
      writer = csv.writer(f, delimiter='\t')
      for taxon_id in taxon_info:
        if taxon_id in present_taxa:
          taxon_percent = (taxon_basepairs[taxon_id] / total_basepairs) * 100
          writer.writerow([
              taxon_percent,
              taxon_basepairs[taxon_id],
              taxon_info[taxon_id]['rank'],
              taxon_id,
              taxon_info[taxon_id]['name']
          ])
    CODE

    # theiacov parsing blocks - percent human, sc2 and target organism
    # capture human percentage
    percentage_human=$(grep "Homo sapiens" ~{samplename}.report_parsed.txt | cut -f 1)
    if [ -z "$percentage_human" ] ; then percentage_human="0" ; fi
    echo "INFO: Percentage human:"
    echo $percentage_human | tee PERCENT_HUMAN

    # capture target org percentage
    if [ ! -z "~{target_organism}" ]; then
      echo "DEBUG: Target org designated: ~{target_organism}"
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
    echo "INFO: Percentage SARS-CoV-2:"
    echo $percentage_sc2 | tee PERCENT_SC2
    echo "INFO: Percentage target organism (~{target_organism}):"
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