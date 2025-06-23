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
    echo "Genus and species: $GENUS $SPECIES"

    # Calculate Assembly Ratio
    if [[ "~{taxon}" == "NA" || -z "$SPECIES" || -z "$GENUS" ]]; then
      echo "Full taxonomy not available" > ASSEMBLY_RATIO
      echo "Full taxonomy not available" > TAXON_ASSEMBLY_RATIO_STDEV
      echo "Full taxonomy not available" > TAXON_GC_ST_DEV
      echo "Full taxonomy not available" > TAXON_GC_MEAN
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
      echo "Assembly standard deviation of taxonomy: ${ref_assem_st_dev}"
      echo "GC standard deviation of taxonomy: ${ref_gc_st_dev}"
      # Metrics Ratios
      assem_ratio=$(python3 -c "print('{:.5f}'.format(~{genome_length} / (float('${assem_mean}') * 1000000)))")
      echo "Calculated assembly ratio: ${assem_ratio}"
      # Calculate Z-score, {observed}-{expected}/{stdev}
      assem_zscore=$(python3 -c "print('{:.5f}'.format((~{genome_length} - (float('${assem_mean}') * 1000000)) / (float('${ref_assem_st_dev}') * 1000000)))")
      # Ratios
      echo "${assem_ratio}" > ASSEMBLY_RATIO
      # StDev
      echo "${ref_assem_st_dev}" > TAXON_ASSEMBLY_RATIO_STDEV
      echo "${ref_gc_st_dev}" > TAXON_GC_ST_DEV
      # Zscore 
      echo "${assem_zscore}" > ASSEMBLY_ZSCORE
      # Mean
      echo "${gc_mean}" > TAXON_GC_MEAN
    else
      echo "Taxon not found in stats file" > ASSEMBLY_RATIO
      echo "Taxon not found in stats file" > TAXON_ASSEMBLY_RATIO_STDEV
      echo "Taxon not found in stats file" > TAXON_GC_ST_DEV
      echo "Taxon not found in stats file" > TAXON_GC_MEAN
      echo "Taxon not found in stats file" > ASSEMBLY_ZSCOREw
    fi
  >>>
  output {
    String read1_raw_q30 = read_string("read1_raw_q30")
    String read2_raw_q30 = read_string("read2_raw_q30")
    String read1_clean_q30 = read_string("read1_clean_q30")
    String read2_clean_q30 = read_string("read2_clean_q30")
    String assembly_ratio = read_string("ASSEMBLY_RATIO")
    String taxon_assembly_ratio_stdev = read_string("TAXON_ASSEMBLY_RATIO_STDEV")
    String taxon_gc_percent_stdev = read_string("TAXON_GC_ST_DEV")
    String taxon_gc_mean = read_string("TAXON_GC_MEAN")
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