version 1.0

task metabuli {
  input {
    File read1
    File? read2
    String samplename
    String? taxon_id
    Int seq_mode = 3 # 1: SE, 2: PE, 3: ONT
    Boolean extract_unclassified = false
    File metabuli_db = "gs://theiagen-public-resources-rp/reference_data/databases/metabuli/refseq_virus-v223.tar.gz"
    File taxdump_path = "gs://theiagen-public-resources-rp/reference_data/databases/metabuli/ncbi_taxdump_20260211.tar.gz"
    Float min_score = 0 # metabuli: Min. sequence similarity score (0.0-1.0) [0.000]
    Float min_sp_score = 0 # metabuli: Min. score for species- or lower-level classification. [0.000]
    Float min_percent_coverage = 0 # metabuli: Min. query coverage (0.0-1.0) [0.000]
    Int cpu = 4
    Int memory = 16
    Int disk_size = 250
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/metabuli:1.1.1"
  }
  command <<<
    set -euo pipefail

    # set status
    metabuli_status=PASS

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

    # Identify available RAM and subtract one for metabuli
    available_ram=$(free -g | sed -nE '2p' | awk '{print $7-1}')
    if [[ $available_ram -le 0 ]]; then
      available_ram=1
    fi
    echo "DEBUG: Available RAM is ${available_ram} GB"

    metabuli classify \
      --seq-mode ~{seq_mode} \
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
      --max-ram ${available_ram}

    # Quantify percent human
    awk '$6 == "Homo" && $7 == "sapiens" {print $1}' output_dir/~{samplename}_report.tsv > PERCENT_HUMAN
    if [[ ! -s PERCENT_HUMAN ]]; then 
      echo "DEBUG: Homo sapiens comprises 0% of reads"
      echo "0" > PERCENT_HUMAN
    else
      echo "DEBUG: Homo sapiens comprises $(cat PERCENT_HUMAN)% of reads"
    fi

    # Extract the reads
    echo "" > PERCENT_TARGET_LINEAGE
    if [[ -n "~{taxon_id}" ]]; then
      echo "DEBUG: Extracting reads"
      awk '$5 == "~{taxon_id}" {print $1}' output_dir/~{samplename}_report.tsv > PERCENT_TARGET_LINEAGE
      if [[ -s PERCENT_TARGET_LINEAGE ]]; then
        echo "DEBUG: Taxon ID ~{taxon_id} found in report, proceeding with read extraction"
        echo "DEBUG: ~{taxon_id} comprises $(cat PERCENT_TARGET_LINEAGE)% of reads"

        metabuli extract \
          ~{read1} ~{read2} \
          ./output_dir/~{samplename}_classifications.tsv \
          ${extracted_db} \
          ~{"--tax-id " + taxon_id} \
          --seq-mode ~{seq_mode} \
          --extract-format 2
    
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
  
      else
        echo "0" > PERCENT_TARGET_LINEAGE
        echo "ERROR: Taxon ID ~{taxon_id} not found in report, skipping read extraction"
        metabuli_status="FAIL; taxon not recovered"
      fi
    fi

    echo $metabuli_status > STATUS
  >>>
  output {
    File metabuli_report = "output_dir/~{samplename}_report.tsv"
    File metabuli_classified = "output_dir/~{samplename}_classifications.tsv"
    File? metabuli_read1_extract = "~{samplename}_~{taxon_id}_extracted_1.fq.gz"
    File? metabuli_read2_extract = "~{samplename}_~{taxon_id}_extracted_2.fq.gz"
    File metabuli_krona_report = "output_dir/~{samplename}_krona.html"
    String metabuli_version = read_string("VERSION")
    String metabuli_status = read_string("STATUS")
    String metabuli_docker = docker
    String metabuli_database = metabuli_db
    String metabuli_percent_target_lineage = read_string("PERCENT_TARGET_LINEAGE")
    Float metabuli_percent_human = read_float("PERCENT_HUMAN")
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