version 1.0

task arln_stats {
  input {
    String samplename
    String taxon
    String workflow_type
    Int genome_length
    File? read1_raw
    File? read2_raw
    File? read1_clean
    File? read2_clean
    Int cpu = 2
    Int memory = 5
    Int disk_size = 10
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/arln_stats:1.0.0"
  }
  command <<<
    set -euo pipefail

    # Generate Q30 statistics raw reads
    if [[ ~{workflow_type} == "pe" ]]; then
      python /scripts/q30.py -i ~{read1_raw} ~{read2_raw} > RAW_Q30
      awk '/read1/ {printf "%.2f\n", $2}' RAW_Q30 > read1_raw_q30
      awk '/read2/ {printf "%.2f\n", $2}' RAW_Q30 > read2_raw_q30
      # Generate Q30 statistics cleaned reads
      python /scripts/q30.py -i ~{read1_clean} ~{read2_clean} > CLEAN_Q30
      awk '/read1/ {printf "%.2f\n", $2}' CLEAN_Q30 > read1_clean_q30
      awk '/read2/ {printf "%.2f\n", $2}' CLEAN_Q30 > read2_clean_q30
    elif [[ ~{workflow_type} == "se" || ~{workflow_type} == "ont" ]]; then
      python /scripts/q30.py -i ~{read1_raw} > RAW_Q30
      awk '/read1/ {printf "%.2f\n", $2}' RAW_Q30 > read1_raw_q30
      touch read2_raw_q30
      python /scripts/q30.py -i ~{read1_clean} > CLEAN_Q30
      awk '/read1/ {printf "%.2f\n", $2}' CLEAN_Q30 > read1_clean_q30
      touch read2_clean_q30
    else
      touch read1_raw_q30
      touch read2_raw_q30
      touch read1_clean_q30
      touch read2_clean_q30
    fi

    GENUS=$(echo ~{taxon} | awk '{print $1}')
    SPECIES=$(echo ~{taxon} | awk '{print $2}')
    echo "$GENUS $SPECIES"

    # Calculate Assembly Ratio
    if [[ "~{taxon}" == "NA" || -z "$SPECIES" ]]; then
      echo "Genus and or Species not found in assembly stats file, check classification level." > assem_ratio_with_stdev
    elif grep -q "~{taxon}" /data/NCBI_Assembly_stats_20240124.txt; then
      assem_mean=$(grep "~{taxon}" /data/NCBI_Assembly_stats_20240124.txt | awk '{print $6}')
      echo "${assem_mean}"
      st_dev=$(grep "~{taxon}" /data/NCBI_Assembly_stats_20240124.txt | awk '{print $7}' | xargs printf "%.5f\n")
      echo "${st_dev}"
      assem_ratio=$(python3 -c "print('{:.4f}'.format((float('${assem_mean}') * 1000000)/ ~{genome_length}))")
      echo "${assem_ratio}x(${st_dev})" > assem_ratio_with_stdev
    else
      echo "Taxon not found in assembly stats file" > assem_ratio_with_stdev
    fi
  >>>
  output {
    String? read1_raw_q30 = read_string("read1_raw_q30")
    String? read2_raw_q30 = read_string("read2_raw_q30")
    String? read1_clean_q30 = read_string("read1_clean_q30")
    String? read2_clean_q30 = read_string("read2_clean_q30")
    String assembly_ratio = read_string("assem_ratio_with_stdev")
    String docker_version = docker
  }
  runtime {
    cpu: cpu
    memory: "~{memory} GB"
    disk: "~{disk_size} GB"
    docker: docker
    maxRetries: 0
  }
}