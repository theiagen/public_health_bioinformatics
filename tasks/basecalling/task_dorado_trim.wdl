version 1.0

task dorado_trim {
  input {
    Array[File] fastq_files  # FASTQ files from the demultiplexing step
    File custom_primers # Custom primers in FASTA format
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.8.3"
    Int cpu = 4
    Int memory = 16
    Int disk_size = 100
  }
  command <<< 
    # Create output directory for trimmed FASTQ files
    output_dir="trimmed_fastqs"
    mkdir -p "$output_dir"

    # Loop through each FASTQ file and apply Dorado trim
    for fq in ~{sep=" " fastq_files}; do

      # Extract the base filename from the path
      filename=$(basename "$fq")
        
      # Perform primer trimming and save output in the designated directory
      dorado trim "$fq" \
        --primer-sequences "~{custom_primers}" \
        --threads ~{cpu} \
        --emit-fastq > "$output_dir/$filename"
    done
  >>>
  output {
    Array[File] trimmed_fastq_files = glob("trimmed_fastqs/*.fastq.gz")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}
