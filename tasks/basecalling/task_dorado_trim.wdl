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
