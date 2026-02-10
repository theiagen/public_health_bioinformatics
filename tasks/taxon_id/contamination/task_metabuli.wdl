version 1.0

task metabuli {
  input {
    File read1
    File? read2
    String samplename
    String? taxon_id
    Boolean extract_unclassified = false
    File metabuli_db = "gs://theiagen-public-resources-rp/reference_data/databases/metabuli/refseq_virus-v223.tar.gz"
    File taxdump_path = "gs://theiagen-public-resources-rp/reference_data/databases/metabuli/new_taxdump.tar.gz"
    Float min_score = 0 # metabuli: Min. sequence similarity score (0.0-1.0) [0.000]
    Float min_sp_score = 0 # metabuli: Min. score for species- or lower-level classification. [0.000]
    Float min_percent_coverage = 0 # metabuli: Min. query coverage (0.0-1.0) [0.000]
    Int cpu = 4
    Int memory = 16
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/metabuli:1.1.1"
  }
  command <<<
    set -euo pipefail

    # get version (there is no --version flag for metabuli)
    echo $(metabuli --help) | awk -F'Version: ' '{print $2}' | awk '{print $1}' | tee VERSION

    # Decompress additional taxonomy files necessary for database search
    echo "DEBUG: Decompressing Metabuli DB"
    mkdir taxdump
    tar -C taxdump/ -xzf ~{taxdump_path}

    # Decompress/extract the database
    mkdir db
    tar -C db/ -xzvf ~{metabuli_db}
    extracted_db=$(ls -d db/*/ | head -n 1)

    # Classify the reads
    echo "DEBUG: Classifying reads"
    metabuli classify \
      --seq-mode 3 \
      ~{read1} \
      ~{read2} \
      ${extracted_db} \
      output_dir \
      "~{samplename}" \
      --taxonomy-path taxdump/ \
      ~{"--min-score " + min_score} \
      ~{"--min-sp-score " + min_sp_score} \
      ~{"--min-cov " + min_percent_coverage} \
      --threads ~{cpu} \
      --max-ram ~{memory}

    # Extract the reads
    if [[ ~{if defined(taxon_id) then "true" else "false"} == "true" ]]; then
      echo "DEBUG: Extracting reads"
      metabuli extract \
        --seq-mode 3 \
        ~{read1} ~{read2} \
        ./output_dir/~{samplename}_classifications.tsv \
        ${extracted_db} \
        ~{"--tax-id " + taxon_id}
  
      # Metabuli extract removes the "gz" suffix, then the final file extension
      read1_basename=$(sed -Ee 's/(.*)\.([^\.]+)$/\1/' <<< $(basename ~{read1} .gz))
  
      # Metabuli extract will create a file in the input directory, which is variable between
      # miniwdl and Terra, so we have to dynamically find it -_-
      find . -type f -name ${read1_basename}_~{taxon_id}.fq -exec mv {} . \;
  
      # Output file name of metabuli extract is static. Compress the extracted reads
      gzip ${read1_basename}_~{taxon_id}.fq
  
      # Update the final extracted reads output file name
      cp ${read1_basename}_~{taxon_id}.fq.gz ~{samplename}_~{taxon_id}_extracted_1.fq.gz
  
      if [[ ~{if defined(read2) then "true" else "false"} == "true" ]]; then
        read2_basename=$(sed -Ee 's/(.*)\.([^\.]+)$/\1/' <<< $(basename ~{read2} .gz))
        find . -type f -name ${read2_basename}_~{taxon_id}.fq -exec mv {} . \;
        gzip ${read2_basename}_~{taxon_id}.fq
        cp ${read2_basename}_~{taxon_id}.fq.gz ~{samplename}_~{taxon_id}_extracted_2.fq.gz
      fi
  
      # Metabuli doesn't have a built-in option for extracting unclassified reads
      if [[ ~{extract_unclassified} == "true" ]]; then
        echo "DEBUG: Extracting unclassified reads"
        grep -P "^0\t" output_dir/~{samplename}_classifications.tsv | cut -f 2 > unclassified_reads.txt
        seqkit grep -f unclassified_reads.txt ~{read1} | gzip > ~{samplename}_unclassified_1.fq.gz
        zcat ${read1_basename}_~{taxon_id}.fq.gz ~{samplename}_unclassified_1.fq.gz | gzip > ~{samplename}_~{taxon_id}_extracted_1.fq.gz
        if [[ ~{if defined(read2) then "true" else "false"} == "true" ]]; then
          seqkit grep -f unclassified_reads.txt ~{read2} | gzip > ~{samplename}_unclassified_2.fq.gz
          zcat ${read2_basename}_~{taxon_id}.fq.gz ~{samplename}_unclassified_2.fq.gz | gzip > ~{samplename}_~{taxon_id}_extracted_2.fq.gz
        fi
    fi
  >>>
  output {
    File metabuli_report = "output_dir/~{samplename}_report.tsv"
    File metabuli_classified = "output_dir/~{samplename}_classifications.tsv"
    File? metabuli_read1_extract = "~{samplename}_~{taxon_id}_extracted_1.fq.gz"
    File? metabuli_read2_extract = "~{samplename}_~{taxon_id}_extracted_2.fq.gz"
    File metabuli_krona_report = "output_dir/~{samplename}_krona.html"
    String metabuli_version = read_string("VERSION")
    String metabuli_docker = docker
    String metabuli_database = metabuli_db
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 0
  }
}