version 1.0

task arln_stats {
    input {
        String samplename
        String taxon
        Int genome_length
        Int cg_bases_r1
        Int cg_bases_r2
        Int cg_bases_cleaned_r1
        Int cg_bases_cleaned_r2
        File cg_report_raw 
        File cg_report_cleaned
        File read1_raw
        File read2_raw
        File read1_clean
        File read2_clean

        Int cpu = 2
        Int memory = 5
        Int disk_size = 10
        String docker = "arln_stats:0.1.0"
    }
    
    command <<<
        set -euo pipefail
        
        # Relay CG Pipeline Basepair Stats
        echo "Processing ARLN statistics"
        echo "~{cg_bases_r1}" | tee cg_bases_r1_out
        echo "~{cg_bases_r2}" | tee cg_bases_r2_out
        echo "~{cg_bases_cleaned_r1}" | tee cg_bases_cleaned_r1_out
        echo "~{cg_bases_cleaned_r2}" | tee cg_bases_cleaned_r2_out

        # Generate Q30 statistics raw reads
        python /scripts/q30.py -i ~{read1_raw} ~{read2_raw} > RAW_Q30
        awk '/read1/ {printf "%.2f\n", $2}' RAW_Q30 > read1_raw_q30
        awk '/read2/ {printf "%.2f\n", $2}' RAW_Q30 > read2_raw_q30
        # Generate Q30 statistics cleaned reads
        python /scripts/q30.py -i ~{read1_clean} ~{read2_clean} > CLEAN_Q30
        awk '/read1/ {printf "%.2f\n", $2}' CLEAN_Q30 > read1_clean_q30
        awk '/read2/ {printf "%.2f\n", $2}' CLEAN_Q30 > read2_clean_q30

        # Calculate Assembly Ratio
        assem_mean=$(grep "~{taxon}" /data/NCBI_Assembly_stats_20240124.txt | awk '{print $6}')
        st_dev=$(grep "~{taxon}" /data/NCBI_Assembly_stats_20240124.txt | awk '{print $7}' | xargs printf "%.4f\n")
        assem_ratio=$(python3 -c "print('{:.4f}'.format((float('${assem_mean}') * 1000000)/ ~{genome_length}))")
        echo "${assem_ratio}x(${st_dev})" > assem_ratio_with_stdev

    >>>
    
    output {
        String cg_bases_r1_out = read_int("cg_bases_r1_out")
        String cg_bases_r2_out = read_int("cg_bases_r2_out")
        String cg_bases_cleaned_r1_out = read_int("cg_bases_cleaned_r1_out")
        String cg_bases_cleaned_r2_out = read_int("cg_bases_cleaned_r2_out")
        Float read1_raw_q30 = read_float("read1_raw_q30")
        Float read2_raw_q30 = read_float("read2_raw_q30")
        Float read1_clean_q30 = read_float("read1_clean_q30")
        Float read2_clean_q30 = read_float("read2_clean_q30")
        String assem_ratio = read_string("assem_ratio_with_stdev")
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