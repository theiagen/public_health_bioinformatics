version 1.0

task dorado_trim {
  input {
    Array[File] fastq_files
    File custom_primers # Custom primers in FASTA format
    
    Int cpu = 4
    Int disk_size = 100
    # this is not the most up-to-date docker image because dorado trim in v0.9.0-cuda12.2.0 is bugged
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.8.3"
    Int memory = 16
  }
  command <<< 
    set -euo pipefail
    
    # Capture Dorado version
    dorado --version 2>&1 | head -n1 | tee DORADO_VERSION

    # Create output directory for trimmed FASTQ files
    output_dir="trimmed_fastqs"
    mkdir -p "$output_dir"

    echo "DEBUG: running dorado trim on the input FASTQ files..."
    # Loop through each FASTQ file and apply Dorado trim
    for fq in ~{sep=" " fastq_files}; do

      # Extract the base filename from the path, and remove .gz extension since dorado trim outputs uncompressed FASTQs
      filename=$(basename "$fq" | sed 's|.gz||')
        
      # Perform primer trimming and save output in the designated directory
      dorado trim "$fq" \
        --primer-sequences "~{custom_primers}" \
        --threads ~{cpu} \
        --emit-fastq > "$output_dir/$filename"
    done
    echo "DEBUG: dorado trim completed on all FASTQ files."

    echo "DEBUG: listing contents of $output_dir:"
    ls -lh "$output_dir"

    # compress FASTQs since dorado trim emits uncompressed FASTQs
    pigz -f  ${output_dir}/*.fastq

    echo "DEBUG: listing contents of $output_dir after pigz compression:"
    ls -lh "$output_dir"
  >>>
  output {
    String dorado_docker = docker
    String dorado_version = read_string("DORADO_VERSION")
    Array[File] trimmed_fastq_files = glob("trimmed_fastqs/*.fastq.gz")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}
