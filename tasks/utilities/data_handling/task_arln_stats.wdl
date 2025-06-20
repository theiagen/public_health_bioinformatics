version 1.0

task arln_stats {
  input {
    String samplename
    String taxon
    String workflow_type
    Int genome_length
    Float gc_percent
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
    if [[ "~{taxon}" == "NA" || -z "$SPECIES" || -z "$GENUS" ]]; then
      echo "Full taxonomy not available" > ASSEMBLY_RATIO
      echo "Full taxonomy not available" > TAXON_ASSEMBLY_RATIO_STDEV
      echo "Full taxonomy not available" > GC_RATIO
      echo "Full taxonomy not available" > TAXON_GC_ST_DEV
      echo "Full taxonomy not available" > GC_ZSCORE
      echo "Full taxonomy not available" > ASSEMBLY_ZSCORE
    elif grep -q "~{taxon}" /data/NCBI_Assembly_stats_20240124.txt; then
      # NCBI Stats Assem Means
      assem_mean=$(grep "~{taxon}" /data/NCBI_Assembly_stats_20240124.txt | awk '{print $6}')
      gc_mean=$(grep "~{taxon}" /data/NCBI_Assembly_stats_20240124.txt | awk '{print $12}')
      echo "Assembly Mean Length: ${assem_mean}"
      echo "GC Mean: ${gc_mean}"
      # NCBI Stats StDevs
      ref_assem_st_dev=$(grep "~{taxon}" /data/NCBI_Assembly_stats_20240124.txt | awk '{print $7}' | xargs printf "%.5f\n")
      ref_gc_st_dev=$(grep "~{taxon}" /data/NCBI_Assembly_stats_20240124.txt | awk '{print $13}' | xargs printf "%.5f\n")
      echo "${ref_assem_st_dev}"
      echo "${ref_gc_st_dev}"
      # Metrics Ratios
      assem_ratio=$(python3 -c "print(~{genome_length} / (float('${assem_mean}') * 1000000))")
      gc_ratio=$(python3 -c "print(~{gc_percent} / float('${gc_mean}'))")
      echo "${assem_ratio}"
      echo "${gc_ratio}"
      # Calculate Z-score, {observed}-{expected}/{stdev}
      assem_zscore=$(python3 -c "print((~{genome_length} - (float('${assem_mean}') * 1000000)) / (float('${ref_assem_st_dev}') * 1000000))")
      gc_zscore=$(python3 -c "print((~{gc_percent} - float('${gc_mean}')) / float('${ref_gc_st_dev}'))")
      # Ratios
      echo "${assem_ratio}" > ASSEMBLY_RATIO
      echo "${gc_ratio}" > GC_RATIO
      # StDev
      echo "${ref_assem_st_dev}" > TAXON_ASSEMBLY_RATIO_STDEV
      echo "${ref_gc_st_dev}" > TAXON_GC_ST_DEV
      # Zscore 
      echo "${assem_zscore}" > ASSEMBLY_ZSCORE
      echo "${gc_zscore}" > GC_ZSCORE 

    else
      echo "Taxon not found in stats file" > ASSEMBLY_RATIO
      echo "Taxon not found in stats file" > TAXON_ASSEMBLY_RATIO_STDEV
      echo "Taxon not found in stats file" > GC_RATIO
      echo "Taxon not found in stats file" > TAXON_GC_ST_DEV
      echo "Taxon not found in stats file" > GC_ZSCORE
      echo "Taxon not found in stats file" > ASSEMBLY_ZSCORE
    fi
  >>>
  output {
    String read1_raw_q30 = read_string("read1_raw_q30")
    String read2_raw_q30 = read_string("read2_raw_q30")
    String read1_clean_q30 = read_string("read1_clean_q30")
    String read2_clean_q30 = read_string("read2_clean_q30")
    String assembly_ratio = read_string("ASSEMBLY_RATIO")
    String taxon_assembly_ratio_stdev = read_string("TAXON_ASSEMBLY_RATIO_STDEV")
    String gc_percent_ratio = read_string("GC_RATIO")
    String taxon_gc_percent_stdev = read_string("TAXON_GC_ST_DEV")
    String gc_zscore = read_string("GC_ZSCORE")
    String assembly_zscore = read_string("ASSEMBLY_ZSCORE")
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